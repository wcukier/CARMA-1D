import numpy as np

import sys
sys.path.insert(0, '..')

from tqdm import tqdm
import carmapy

from matplotlib import pyplot as plt


Teff_grid = range(900, 2401, 100)
g_grid = [31, 100, 316, 1000, 3160]
# g_grid = [316]
# m_grid = [0.0]
m_grid = [-0.5, 0.0, 0.5]

numdens = np.zeros((16, 5,3, 90, 11, 80, 99))

for i_t in tqdm(range(len(Teff_grid))):
    for i_g in range(len(g_grid)):
        for i_m in range(len(m_grid)):

            teff = Teff_grid[i_t]
            g = g_grid[i_g]
            m = m_grid[i_m]
            
            # for _ in range(1):
            try:
                carma = carmapy.load_carma(f"../outputs/t{teff}g{g}m{m:+.1f}")
                carma.read_results()
                # if not numdems:
                #     numdens = np.zeros((len(Teff_grid), len(g_grid), len(m_grid), *carma.results.numden.shape)) * np.nan
                n_t = carma.results.numden.shape[3]
                numdens[i_t, i_g, i_m, :, :, :, :n_t] = carma.results.numden
            except Exception as e:
                print(f"{i_t}, {i_g}, {i_m}, :{e}")
    
np.save("numdens.npy", numdens)