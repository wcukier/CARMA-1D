import numpy as np
from .constants import *

def load_results(carma):
    path = carma.name
    path_end = path.split("/")[-1]
    file_path = path+f"/bd {path_end}.txt"
    
    f = open(file_path)
    NZ, NGROUP, NELEM, NBIN, NGAS, nstep, iskip = np.array(f.readline().split(),
                                                           dtype=int)

    
    if ((NZ != carma.NZ) + 
        (NGROUP != len(carma.groups))+
        (NELEM != len(carma.elems))+
        (NBIN != carma.NBIN) +
        (NGAS != len(carma.gasses))+
        (nstep - 1 != carma.n_tstep)+
        (iskip != carma.output_gap)
    ):
        raise ValueError(f"Output file inconsistent with carma run")
    
    r = np.zeros(NBIN)
    rmass = np.zeros((NBIN, NGROUP))
    
    for i in range(NGROUP):
        for j in range(NBIN):
            _, _, r[j], rmass[j, i], _, _, _ = np.array(f.readline().split(), dtype=float)
        
    
    kzz = np.zeros(NZ)
    P = np.zeros(NZ)
    T = np.zeros(NZ)
    Z = np.zeros(NZ)

    for i in range(NZ):
        _, Z[i], _, P[i], T[i], kzz[i] = np.array(f.readline().split(), 
                                                  dtype=float)

    f.readline()
    f.readline()

    for j in range(NBIN):
        for i in range(NZ):
            f.readline()
        
    NT = int(nstep/iskip)

    numden = np.zeros((NZ, NELEM, NBIN, NT))
    gas_abund = np.zeros((NZ, NGAS, NT))
    sat_vp = np.zeros((NZ, NGAS, NT))
    ts = np.zeros(NT)
    
    for it in range(NT):
        t_step = f.readline()
        if t_step:
            ts[it] = t_step
            for ibin in range(NBIN):
                for iz in range(NZ):
                    line = np.array(f.readline().split(), dtype=float)
                    for ielem in range(NELEM):
                        numden[iz, ielem, ibin, it] = line[ielem+2]
                    for igas in range(NGAS):
                        gas_abund[iz, igas, it] = line[NELEM + 2+ 2*igas]
                        sat_vp[iz, igas, it] = line[NELEM + 3 + 2*igas]
        else:
            break   
        
    numden_groups = np.zeros((NZ, NGROUP, NBIN, NT))
    
    
    for i, key in enumerate(carma.groups.keys()):
        group = carma.groups[key]
        if group.mantle:
            numden_groups[:, group.igroup-1, :, :] = numden[:, group.mantle.ielem-1, :, :]
        else:
            numden_groups[:, group.igroup-1, :, :] = numden[:, group.core.ielem-1, :, :]
        if np.any( numden_groups[:, group.igroup-1, :, :] < 0):
            print(key, group.mantle, np.max( numden[:, group.core.ielem-1, :, :]))
            raise
    
        
    carma.results['rmass'] = rmass
    carma.results['r'] = r
    carma.results["numden"] = numden_groups[:,:,:,:it]
    carma.results["gas_abund"] = gas_abund[:,:,:it]
    carma.results["sat_vp"] = sat_vp[:,:,:it]
    carma.results["ts"] = ts[:it]
                    

