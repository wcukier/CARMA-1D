import carmapy
import numpy as np
from scipy.interpolate import LinearNDInterpolator
import pyfastchem
import sys

k_B = 1.381e-16 #erg/K
BAR_TO_BARYE = 1e6
m_p = 1.673e-24 # g
STEFAN_CONSTANT = 5.67e-5 #erg/(cm^2 * s * K^4)
N_A = 6.022e23

atomic_mass = np.array([0, 2, 1, 1, 1, 66.939, 63.866, 44.009, 4.002, 18.015, 16.04, 28.01, 17.031, 28.014, 33.998, 34.08, 55.845, 22.990, 39.098])


g = 316 #cm/s^2
# teff = 2400 #K
f = 4
m = "0.0"

def get_fastchem_abundances(T, P, species):

  temperature = T
  pressure = np.array(P) / BAR_TO_BARYE

  fastchem = pyfastchem.FastChem(
    'inputs/fastchem/asplund_2009_extended.dat',
    'inputs/fastchem/logK.dat',
    1)


  input_data = pyfastchem.FastChemInput()
  output_data = pyfastchem.FastChemOutput()

  input_data.temperature = temperature
  input_data.pressure = pressure
  fastchem_flag = fastchem.calcDensities(input_data, output_data)

  print("FastChem reports:")
  print("  -", pyfastchem.FASTCHEM_MSG[fastchem_flag])

  if np.amin(output_data.element_conserved[:]) == 1:
    print("  - element conservation: ok")
  else:
    print("  - element conservation: fail")
    
  number_densities = np.array(output_data.number_densities)
  nmr = number_densities / np.repeat((P/(k_B * T))[:, np.newaxis], number_densities.shape[1], axis=1)
  
  print(number_densities.shape)
  ret = []
  for s in species:
    index = fastchem.getGasSpeciesIndex(s)
    if index == pyfastchem.FASTCHEM_UNKNOWN_SPECIES:
      raise ValueError(f"{s} is an unknown species")
    ret.append(nmr[:, index])
  
  
  return np.array(ret)



def load_profiles(teff, g, f, m):
    run_num = sys.argv[1]
    
    
    
    data = np.genfromtxt(f"inputs/pressure-temperature_profiles/t{teff}g{g}f{f}_m{m}_co1.0.pt", skip_header=2)
    data[:, 1] = data[:, 1] * BAR_TO_BARYE

    P_levels = data[:, 1]
    T_levels = data[:, 2]

    P_centers = (P_levels[:-1] + P_levels[:1])/2
    T_centers = (T_levels[:-1] + T_levels[:1])/2



    μ_levels = np.sum(atomic_mass * data[:, 3:], axis = 1)/np.sum(data[:, 3:], axis = 1)
    μ_centers = 1/2 * (μ_levels[1:] + μ_levels[:-1])

    H_levels = k_B * T_levels/(μ_levels * m_p * g)
    H_centers = k_B * T_centers/(μ_centers * m_p * g)

    z_levels = np.zeros(data.shape[0])

    for i in range(1, z_levels.shape[0]):
        delta_z = H_levels[i] * np.log(data[i, 1]/data[i-1, 1])

        z_levels[i] = z_levels[i-1] + delta_z

    z_centers = (z_levels[1:] + z_levels[:-1])/2


    c_p = 5/2 * k_B/(μ_centers * m_p) # erg / (g K)
    ρ = μ_centers * P_centers * m_p/(k_B * T_centers) 
    F = STEFAN_CONSTANT * teff**4

    Γ = (T_levels[1:] - T_levels[:-1]) / (z_levels[1:] - z_levels[:-1])
    Γ_a = g/c_p
    lapse_ratio = Γ/Γ_a

    dlnP = np.log(P_levels[:-1]/P_levels[1:])
    lapse_ratio = (T_levels[:-1]- T_levels[1:])/dlnP/(2/7 * T_centers)


    L = H_centers * np.array([max(0.1, lapse_ratio[i]) for i in range(len(Γ))])
    k_zz = H_centers/3 * (L/H_centers)**(4/3)* (k_B * F/(μ_centers*m_p * ρ * c_p))**(1/3)

    k_zz_levels = np.zeros(P_levels.shape)
    k_zz_levels[1:-1] = (k_zz[1:] + k_zz[:-1])/2
    k_zz_levels[0] = k_zz_levels[1] + 2* (k_zz[0] - k_zz_levels[1])
    k_zz_levels[-1] = k_zz_levels[-2] + 2* (k_zz[-1] - k_zz_levels[-2])


    z_centers = np.flip(z_levels[-1] - z_centers)
    z_levels =np.flip(z_levels[-1] - z_levels)
    P_centers = np.flip(P_centers)
    P_levels = np.flip(P_levels)
    T_centers = np.flip(T_centers)
    k_zz_levels = np.flip(k_zz_levels)
    return z_centers, z_levels, P_centers, P_levels, T_centers, k_zz_levels
    

if __name__ == "__main__":
    run_num = int(sys.argv[1])
    
    Teff_grid = range(1400, 2401, 100)
    teff = Teff_grid[run_num]
    print(f"Teff: {teff}")
    
    (z_centers, z_levels, P_centers, 
     P_levels, T_centers, k_zz_levels) = load_profiles(teff, g, f, m)
    
    species = ["Ti", "Fe", "Mg", "Mg", "Cr", "Mn", "Na", "Zn", "Al"]
    abund = get_fastchem_abundances(T_centers, P_centers, species)
    
    carma = carmapy.Carma(f"outputs/{teff}")
    carma.add_gas("H2O")


    carma.set_stepping(dt=100, output_gap=1000, n_tstep=100000)
    carma.set_physical_params(surface_grav=316,
                            wt_mol = 2.389,
                            r_planet = 6.991e9)

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

    nmr_dict = {
        "H2O": 0,
        "TiO2": abund[0, 0],
        "Fe": abund[1, 0],
        "Mg2SiO4": abund[2, 0],
        "Cr": abund[3, 0],
        "MnS": abund[4, 0],
        "Na2S": abund[5, 0], 
        "ZnS": abund[6, 0], 
        "KCl": abund[7, 0], 
        "Al2O3": abund[8, 0]
    }

    carma.set_nmr(nmr_dict)

    carma.add_z(z_centers, z_levels)
    carma.add_P(P_centers, P_levels)
    carma.add_kzz(k_zz_levels)
    carma.add_T(T_centers)
    
    carma.run()