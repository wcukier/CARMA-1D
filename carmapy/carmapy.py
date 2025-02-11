from .constants import *
import os
import shutil
import f90nml
import numpy as np
import subprocess
import warnings

class Carma:
    def __init__(self, name):
        
        self.NZ = 0
        self.P_levels = None
        self.P_centers = None
        self.z_levels = None
        self.z_centers = None
        self.T_centers = None
        self.kzz_levels = None
        

        self.groups = {}
        self.growth = []
        self.elems = {}
        self.gasses = {}
        self.nucs = []
        self.name = name
        self.surface_grav = 3160
        self.wt_mol = 2.3
        self.r_planet = 6.991e9
        self.restart = False
        
        self.dt = 1000
        self.output_gap = 1000
        self.n_tstep = 1_000_000

        
    def set_stepping(self, dt=None, output_gap = None, n_tstep = None):
        if dt:
            if dt != int(dt):
                raise TypeError("dt must be a integer")
            dt = int(dt)
            if dt < 0:
                raise ValueError("dt must be positive")
        
        if output_gap:
            if output_gap != int(output_gap):
                raise TypeError("dt must be a integer")
            output_gap = int(output_gap)
            if dt <= 0:
                raise ValueError("output_gap must be positive")
            
        if n_tstep:
            if n_tstep != int(n_tstep):
                raise TypeError("n_tstep must be a integer")
            n_tstep = int(n_tstep)
            if n_tstep < 0:
                raise ValueError("n_tstep must be positive")
            
        if dt: self.dt = dt
        if output_gap: self.output_gap = output_gap
        if n_tstep: self.n_tstep = n_tstep
            
        
    def set_physical_params(self, surface_grav = None, wt_mol = None, r_planet = None, use_jovian_radius=False):
        if surface_grav:
            if surface_grav < 0:
                raise ValueError("Surface Gravity must be positive")
       
            
        if wt_mol:
            if wt_mol < 0:
                raise ValueError("Molar Weight must be positive")
            if wt_mol > 3 or wt_mol < 2:
                warnings.warn(f"Typical values of wt_mol are between 2 and 3.  Your value is {wt_mol}.")
           
            
        if r_planet:
            if r_planet < 0:
                raise ValueError("Planet Radius must be positive")
            if r_planet < 20 and not use_jovian_radius:
                warnings.warn("You specified a planetary radius under 20.  Assuming you meant in units of Jovian radius.  \n set use_jovian_radius=True to supress this warning")
            if r_planet > 1e3 and use_jovian_radius:
                raise ValueError(f"The specified planetary radius of {r_planet} jovian radii is too high")
            
        
        if surface_grav: self.surface_grav = surface_grav
        if wt_mol: self.wt_mol = wt_mol
        if r_planet:
            if use_jovian_radius:
                self.r_planet = r_planet * JUPITER_RADIUS
            else:
                self.r_planet = r_planet
   
    
    def add_kzz(self,levels):
        if self.NZ:
            if len(levels) != self.NZ + 1:
                raise ValueError(f"centers must be {self.NZ+1} long to be compatible with other input.\nYour data was {len(levels)} long.")
        else:
            self.NZ = len(levels) - 1
        self.kzz_levels = levels
    
    def add_P(self, centers, levels):
        if self.NZ:
            if len(centers) != self.NZ:
                raise ValueError(f"centers must be {self.NZ} long to be compatible with other input.\nYour data was {len(centers)} long.")
            if len(levels) != self.NZ + 1:
                raise ValueError(f"centers must be {self.NZ+1} long to be compatible with other input.\nYour data was {len(levels)} long.")
        else:
            if len(levels) != len(centers) + 1:
                raise ValueError(f"levels must be one longer than centers.  Your centers was {len(centers)} and your levels was {len(levels)}")
            self.NZ = len(centers)
        self.P_centers = centers
        self.P_levels = levels
        
    def add_T(self, centers):
        if self.NZ:
            if len(centers) != self.NZ:
                raise ValueError(f"centers must be {self.NZ} long to be compatible with other input.\nYour data was {len(centers)} long.")
        else:
           self.NZ = len(centers)
        self.T_centers = centers

        
    def add_z(self, centers, levels):
        if self.NZ:
            if len(centers) != self.NZ:
                raise ValueError(f"centers must be {self.NZ} long to be compatible with other input.\nYour data was {len(centers)} long.")
            if len(levels) != self.NZ + 1:
                raise ValueError(f"centers must be {self.NZ+1} long to be compatible with other input.\nYour data was {len(levels)} long.")
        else:
            if len(levels) != len(centers) + 1:
                raise ValueError(f"levels must be one longer than centers.  Your centers was {len(centers)} and your levels was {len(levels)}")
            self.NZ = len(centers)
        self.z_centers = centers
        self.z_levels = levels

    
    def add_het_group(self, gas,  seed_group, rmin, mucos=None):
        if type(gas) == type(""):
            gas = self.gasses.get(gas, Gas(gas, len(self.gasses) + 1))
            self.gasses[gas.name] = gas
        
        if type(seed_group) == type(""):
            seed_group = self.groups["Pure "+seed_group.split(" ")[-1]]
        
        name =  gas.name + " on " + seed_group.name.split(" ")[-1]
        group = Group(len(self.groups)+1, name, rmin)
        self.groups[name] = group
   
        


        if not mucos:
            mucos = mucos_dict[gas.name][seed_group.name.split(" ")[-1]]
        
        self.nucs.append(Nuc(seed_group, group, True, gas, mucos))




        
        mantle_elem = Element(gas.name + " Mantle", len(self.elems)+1, 
                              group, cond_rho[gas.name], "Volatile", 
                              igas_dict[gas.name])
        self.elems[mantle_elem.name] = mantle_elem
        group.mantle = mantle_elem

        core_elem = seed_group.coreify(len(self.elems)+1, group, gas.name)
        self.elems[core_elem.name] = core_elem
        group.core = core_elem
        
        growth = Growth(mantle_elem, gas)
        self.growth.append(growth)
        
    def add_hom_group(self, gas, rmin):
        if type(gas) == type(""):
            gas = self.gasses.get(gas, Gas(gas, len(self.gasses) + 1))
            self.gasses[gas.name] = gas
            
        name = "Pure "+ gas.name
        group = Group(len(self.groups)+1, name, rmin)
        self.groups[name] = group
        

        
        self.nucs.append(Nuc(group, group, False, gas,  0))
        
        elem = Element("Pure "+ gas.name, len(self.elems)+1, 
                            group, cond_rho[gas.name], "Volatile", 
                            igas_dict[gas.name])
        group.core = elem
        self.elems[elem.name] = elem
        self.growth.append(Growth(elem, gas))
    
    def add_gas(self, gas, **kwargs):
        self.gasses[gas] = self.gasses.get(gas, Gas(gas, len(self.gasses)+1, **kwargs))
      
    def set_nmr(self, nmr_dict):
        for key in nmr_dict.keys():
            self.gasses[key].nmr = nmr_dict[key]
    
    def run(self, path=None):
        if not path: path = self.name
        
        os.makedirs(path, exist_ok=True)
        os.makedirs(path+"/inputs", exist_ok=True)
        shutil.copy("carmapy/diamondback_test.exe", path)
        
        nml = {
            "io_files": {
                "filename": self.name,
                "filename_restart": self.name+"_restart",
                "fileprefix": "bd",
                "gas_input_file": "inputs/gas_input.txt",
                "centers_file": "inputs/centers.txt",
                "levels_file": "inputs/levels.txt",
                "groups_file": "inputs/groups.txt",
                "elements_file": "inputs/elements.txt",
                "gases_file": "inputs/gasses.txt",
                "growth_file": "inputs/growth.txt",
                "nuc_file": "inputs/nucleation.txt"
                },
            "physical_params" : {
                "wtmol_air_set": self.wt_mol,
                "grav_set": self.surface_grav,
                "rplanet":  self.r_planet
                },
            "input_params": {
                "NZ": self.NZ-1,
                "NELEM": len(self.elems),
                "NGROUP": len(self.groups),
                "NGAS": len(self.gasses),
                "NBIN": 80,
                "NSOLUTE": 1,
                "NWAVE": 0,
                "ilongitude": 64,
                "irestart": int(self.restart),
                "idiag": 0,
                "iskip": self.output_gap,
                "nstep": self.n_tstep,
                "dtime": self.dt,
                "NGROWTH": len(self.growth),
                "NNUC": len(self.nucs)
                }        
            }
        nml = f90nml.Namelist(nml)
        nml.write(path+"/inputs/input.nml", force=True)
        
        with open(path+"/inputs/groups.txt", "w+") as f:
            for key in self.groups.keys():
                name = '"'+key + '"'
                f.write(f'{name:24s}{self.groups[key].rmin:.15e}\n')
        
        with open(path+"/inputs/gasses.txt", "w+") as f:
            for key in self.gasses.keys():
                name = '"'+key + ' Vapor"'
                f.write(f'{name:24s}{self.gasses[key].wtmol:<.6e}\t{self.gasses[key].ivaprtn:2d}\t{self.gasses[key].icomp:2d}\t{self.gasses[key].wtmol_dif:3.4f}\n')
        
        with open(path+"/inputs/elements.txt", "w+") as f:
            for key in self.elems.keys():
                name = '"'+key + '"'
                proc = '"'+self.elems[key].proc + '"'
                f.write(f'{self.elems[key].group.igroup}\t{name:24s}{self.elems[key].rho:2.4f}\t{proc:15s}\t{self.elems[key].igas:2d}\n')
        
        
        
        with open(path+"/inputs/nucleation.txt", "w+") as f:
            for nuc in self.nucs:
                igas = nuc.gas.igas
                if nuc.is_het:
                    ele_from = nuc.group_from.core.ielem
                    ele_to = nuc.group_to.core.ielem
                    f.write(f'{ele_from:3d}\t{ele_to:3d}\t1\t{igas:3d}\t{ele_from:3d}\t{nuc.mucos:1.8f}\n')

                else:
                    ele_from = nuc.group_from.core.ielem
                    ele_to = ele_from
                    f.write(f'{ele_from:3d}\t{ele_to:3d}\t0\t{igas:3d}\t{0:3d}\t{0:1.8f}\n')
        
        with open(path+"/inputs/growth.txt", "w+") as f:
            for g in self.growth:
                f.write(f"{g.elem.ielem}\t {g.gas.igas}\n")

        with open(path+"/inputs/centers.txt", "w+") as f:
            for i in range(self.NZ):
                f.write(f"{self.z_centers[i]/100}\t{self.P_centers[i]/10}\t{self.T_centers[i]}\n")
        
        with open(path+"/inputs/levels.txt", "w+") as f:
            for i in range(self.NZ):
                f.write(f"{self.z_levels[i]/100}\t{self.P_levels[i]/10}\t{self.kzz_levels[i]}\n")
        
        with open(path+"/inputs/gas_input.txt", "w+") as f:
            for key in self.gasses.keys():
                g = self.gasses[key]
                if g.nmr < 0:
                    raise AttributeError(f"The nmr for {g.name} was not set.")
                f.write(f"{g.nmr:10e}\t")
            f.write("\n")
            for i in range(1, self.NZ):
                for key in self.gasses.keys():
                    g = self.gasses[key]
                    if len(np.shape(g.nmr)) > 0:
                        if len(g.nmr) != self.NZ:
                            raise ValueError(f"The array for nmr of {g.name} is {len(g.nmr)}.  It should be {self.NZ}.")
                        f.write(f"{g.nmr:10e}\t")
                    else:
                        f.write(f"{0.:10e}\t")
                f.write("\n")
        
        wd = os.getcwd()
        os.chdir(path)
        try:
            subprocess.run(["export", "OMP_NUM_THREADS=1"], shell=True,stdout=subprocess.PIPE)
            subprocess.run(["export", "KMP_STACKSIZE=128M"], shell=True,stdout=subprocess.PIPE)
            p = subprocess.Popen("./diamondback_test.exe", shell=True, stdout=subprocess.PIPE)
            
            while p.poll() is None:
                l = p.stdout.readline() # This blocks until it receives a newline.
                print(l.decode('UTF-8'))
            # When the subprocess terminates there might be unconsumed output 
            # that still needs to be processed.
            print(p.stdout.read().decode('UTF-8'))
        except Exception as e:
            print(e)
        os.chdir(wd)
            
            
        
    
    
class Element:
    def __init__(self, name, ielem, group, rho, proc, igas):
        self.name = name
        self.ielem = ielem
        self.group = group
        self.rho = rho
        self.proc = proc
        self.igas = igas
    
class Gas:
    def __init__(self, gas_name, igas, **kwargs):
        self.name = gas_name
        self.igas = igas
        self.wtmol = kwargs.get("wtmol", wtmol_dict[gas_name])
        self.ivaprtn = kwargs.get("ivaprtn",vaprtn_dict[gas_name])
        self.icomp = gcomp_dict[gas_name]
        self.wtmol_dif = kwargs.get("wtmol_dif", wtmol_dif_dict[gas_name])
        self.nmr = kwargs.get("nmr", -1)
        
    
class Nuc:
    def __init__(self, group_from, group_to, is_het, gas, mucos):
        self.group_from = group_from
        self.group_to = group_to
        self.is_het = is_het
        self.gas = gas
        self.ievp2elem = group_from.core
        self.mucos = mucos
    
class Growth:
    def __init__(self, elem, gas):
        self.elem = elem
        self.gas = gas
    
class Group:
    def __init__(self, igroup, name, rmin):
        self.igroup = igroup
        self.name = name
        self.rmin = rmin
        self.core = None
        self.mantle = None
    
    def coreify(self, ielem, group, gas_name=""):
        core_elem = self.core
        
        name = core_elem.name
        name = name.split(" ")[-1]
        name = name + f" Core ({gas_name})"
        
        elem = Element(name, ielem, group, core_elem.rho, "Core Mass", core_elem.igas)
        return elem

def default_carma(name):
    carma = Carma(name)
    carma.add_gas("H2O")

    # Optional, here to preserve ordering
    carma.add_gas("TiO2")
    carma.add_gas("Fe")
    carma.add_gas("Mg2SiO4")
    carma.add_gas("Cr")
    carma.add_gas("MnS")
    carma.add_gas("Na2S")
    carma.add_gas("ZnS")
    carma.add_gas("KCl")
    carma.add_gas("Al2O3")



    carma.add_hom_group("TiO2", 1e-8)
    carma.add_het_group("Al2O3", "TiO2", 1e-8 * 2**(1/3))
    carma.add_het_group("Fe", "TiO2", 1e-8 * 2**(1/3))
    carma.add_het_group("Mg2SiO4", "TiO2", 1e-8 * 2**(1/3))
    carma.add_het_group("Cr", "TiO2", 1e-8 * 2**(1/3))
    carma.add_het_group("MnS", "TiO2", 1e-8 * 2**(1/3))
    carma.add_het_group("Na2S", "TiO2", 1e-8 * 2**(2/3))
    carma.add_hom_group("Fe", 1e-8)
    carma.add_hom_group("Cr", 1e-8)
    carma.add_hom_group("KCl", 1e-8)
    carma.add_het_group("ZnS", "KCl", 1e-8 * 2**(1/3))
    return carma