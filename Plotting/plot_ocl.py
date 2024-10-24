#! /usr/bin/env python2.7
import sys,matplotlib
from matplotlib import colors, ticker, cm
from matplotlib.colors import LogNorm
sys.path.append('../../utils/')
matplotlib.use('Agg')
from numpy import *
from pylab import *
import getopt
import os, re
import pdb
import math as math
import matplotlib.gridspec as gridspec
import colormaps as cmaps


plt.rcParams['figure.figsize'] = (6,5)
plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['legend.borderpad'] = 0.1
plt.rcParams['legend.labelspacing'] = 0.1
plt.rcParams['legend.handletextpad'] = 0.1
plt.rcParams['font.family'] = 'stixgeneral'
plt.rcParams['font.size'] = 18
plt.rcParams['xtick.major.size'] = 3.5
plt.rcParams['xtick.minor.size'] = 2
plt.rcParams['xtick.major.width'] = 0.8
plt.rcParams['xtick.minor.width'] = 0.6
plt.rcParams['ytick.major.size'] = 3.5
plt.rcParams['ytick.minor.size'] = 2
plt.rcParams['ytick.major.width'] = 0.8
plt.rcParams['xtick.minor.width'] = 0.6

wave, ocl_10 = np.loadtxt('nadir_spectra/ocl_1000.txt', unpack = True)
wave, ocl_11 = np.loadtxt('nadir_spectra/ocl_1100.txt', unpack = True)
wave, ocl_12 = np.loadtxt('nadir_spectra/ocl_1200.txt', unpack = True)
wave, ocl_13 = np.loadtxt('nadir_spectra/ocl_1300.txt', unpack = True)
wave, ocl_14 = np.loadtxt('nadir_spectra/ocl_1400.txt', unpack = True)
wave, ocl_15 = np.loadtxt('nadir_spectra/ocl_1500.txt', unpack = True)
wave, ocl_16 = np.loadtxt('nadir_spectra/ocl_1600.txt', unpack = True)
wave, ocl_17 = np.loadtxt('nadir_spectra/ocl_1700.txt', unpack = True)
wave, ocl_18 = np.loadtxt('nadir_spectra/ocl_1800.txt', unpack = True)
wave, ocl_19 = np.loadtxt('nadir_spectra/ocl_1900.txt', unpack = True)
wave, ocl_20 = np.loadtxt('nadir_spectra/ocl_2000.txt', unpack = True)
wave, ocl_21 = np.loadtxt('nadir_spectra/ocl_2100.txt', unpack = True)
wave, ocl_22 = np.loadtxt('nadir_spectra/ocl_2200.txt', unpack = True)
wave, ocl_23 = np.loadtxt('nadir_spectra/ocl_2300.txt', unpack = True)
wave, ocl_24 = np.loadtxt('nadir_spectra/ocl_2400.txt', unpack = True)

z, p_10, t_10 = np.genfromtxt('../run/carma/t1000g1000nc_m0.0.atm', skip_header = 1, skip_footer = 93, unpack = True)


plt.loglog(wave,ocl_10, lw=2, c='r', alpha = .3)
plt.loglog(wave,ocl_11, lw=2, c='r', alpha = .3)
plt.loglog(wave,ocl_12, lw=2, c='r', alpha = .3)
plt.loglog(wave,ocl_13, lw=2, c='r', alpha = .3)
plt.loglog(wave,ocl_14, lw=2, c='r', alpha = .3)
plt.loglog(wave,ocl_15, lw=2, c='r', alpha = .3)
plt.loglog(wave,ocl_16, lw=2, c='r', alpha = .3)
plt.loglog(wave,ocl_17, lw=2, c='r', alpha = .3)
plt.loglog(wave,ocl_18, lw=2, c='r', alpha = .3)
plt.loglog(wave,ocl_19, lw=2, c='r', alpha = .3)
plt.loglog(wave,ocl_20, lw=2, c='r', alpha = .3)
plt.loglog(wave,ocl_21, lw=2, c='r', alpha = .3)
plt.loglog(wave,ocl_22, lw=2, c='r', alpha = .3)
plt.loglog(wave,ocl_23, lw=2, c='r', alpha = .3)
plt.loglog(wave,ocl_24, lw=2, c='r', alpha = .3)
plt.xlim(min(wave), max(wave))
plt.ylim(max(p_10*1e-6), min(p_10*1e-6))
plt.xlabel(r'Wavelength [$\mu$m]')
plt.ylabel(r'Pressure [bar]')
plt.savefig('ocl.png', format = 'png', bbox_inches='tight')
plt.close()

idx = (np.abs(wave - 10)).argmin()

plot1 = [ocl_10[0], ocl_11[0], ocl_12[0], ocl_13[0], ocl_14[0], ocl_15[0], ocl_16[0], ocl_17[0], ocl_18[0], ocl_19[0], ocl_20[0], ocl_21[0], ocl_22[0], ocl_23[0], ocl_24[0]]
plot2 = [ocl_10[len(wave)/2], ocl_11[len(wave)/2], ocl_12[len(wave)/2], ocl_13[len(wave)/2], ocl_14[len(wave)/2], ocl_15[len(wave)/2], ocl_16[len(wave)/2], ocl_17[len(wave)/2], ocl_18[len(wave)/2], ocl_19[len(wave)/2], ocl_20[len(wave)/2], ocl_21[len(wave)/2], ocl_22[len(wave)/2], ocl_23[len(wave)/2], ocl_24[len(wave)/2]]
plot3 = [ocl_10[idx], ocl_11[idx], ocl_12[idx], ocl_13[idx], ocl_14[idx], ocl_15[idx], ocl_16[idx], ocl_17[idx], ocl_18[idx], ocl_19[idx], ocl_20[idx], ocl_21[idx], ocl_22[idx], ocl_23[idx], ocl_24[idx]]

temp = np.linspace(1000, 2400, num=15)

plt.scatter(temp,plot1,s=35, c='#66c2a5')
plt.plot(temp,plot1, lw=2, c='#66c2a5', label = str(round(wave[0],2))+' $\mu$m')
plt.scatter(temp,plot2,s=35, c='#fc8d62')
plt.plot(temp,plot2, lw=2, c='#fc8d62', label = str(round(wave[len(wave)/2],2))+' $\mu$m')
plt.scatter(temp,plot3,s=35, c='#8da0cb')
plt.plot(temp,plot3, lw=2, c='#8da0cb', label = str(round(wave[idx],2))+' $\mu$m')

plt.yscale('log')
plt.ylim(max(p_10*1e-6), min(p_10*1e-6))
plt.xlim(900,2500)
plt.xlabel(r'Equilibrium Temperature [K]')
plt.ylabel(r'Pressure [bar]')
plt.legend(loc='best')
plt.savefig('ocl_temp_trend.png', format = 'png', bbox_inches='tight')
plt.close()

