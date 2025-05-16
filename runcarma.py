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


# g = 316 #cm/s^2
# teff = 2400 #K
f = 4
# m = "0.0"


def load_profiles(teff, g, f, m):
    m_string = f"{m:.1f}"
    if m > 0: m_string = "+" + m_string
    
    
    data = np.genfromtxt(f"inputs/pressure-temperature_profiles/t{teff}g{int(g/100)}f{f}_m{m_string}_co1.0.pt", skip_header=2)

    P_levels = data[:, 1] * BAR_TO_BARYE
    T_levels = data[:, 2]

    μ_levels = np.sum(atomic_mass * data[:, 3:], axis = 1)/np.sum(data[:, 3:], axis = 1)

    k_zz = np.genfromtxt(f"inputs/kzz/t{teff}g{int(g/100)}f{f}_m{m_string}_co1.0.txt")
    k_zz_levels = np.zeros(P_levels.shape)
    k_zz_levels[1:-1] = (k_zz[1:] + k_zz[:-1])/2
    k_zz_levels[0] = k_zz_levels[1] + 2* (k_zz[0] - k_zz_levels[1])
    k_zz_levels[-1] = k_zz_levels[-2] + 2* (k_zz[-1] - k_zz_levels[-2])

    P_levels = np.flip(P_levels)
    T_levels = np.flip(T_levels)
    k_zz_levels = np.flip(k_zz_levels)
    μ_levels = np.flip(μ_levels)

    return P_levels, T_levels, k_zz_levels, μ_levels
    

if __name__ == "__main__":
    run_num = int(sys.argv[1])
    
    Teff_grid = range(900, 2401, 100)
    g_grid = [3100, 10000, 31600, 100000, 316000]
    m_grid = [-0.5, 0.0, 0.5]
    
    
    teff = Teff_grid[run_num % len(Teff_grid)]
    temp_run_num = int(run_num/len(Teff_grid))
    
    g = g_grid[temp_run_num % len(g_grid)]
    temp_run_num = int(temp_run_num/len(g_grid))
    
    m = m_grid[temp_run_num % len(m_grid)]
    
    
    print(f"Teff: {teff}, g:{g}, m:{m}", file=sys.stderr)
    
    P_levels, T_levels, k_zz_levels, μ_levels = load_profiles(teff, g, f, m)
    
    
    carma = carmapy.Carma(f"/scratch/midway3/wcukier/outputs-r4/t{teff}g{int(g/100)}m{m:+.1f}")
    carma.add_gas("H2O")

    carma.set_stepping(dt=10, output_gap=10000, n_tstep=1000000)
    carma.set_physical_params(surface_grav=g,
                            wt_mol = μ_levels[0],
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

    carma.add_P(P_levels)
    carma.add_T(T_levels)
    carma.add_kzz(k_zz_levels)
    
    carma.calculate_z(μ_levels)
    carma.extend_atmosphere(1e10)
    carmapy.populate_fastchem_abundances(carma, metalicity=10**m)

    carma.dt = 100
    carma.output_gap = 1000
    carma.n_tstep = 100_000

    carma.run()