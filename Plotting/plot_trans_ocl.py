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
from scipy.interpolate import interp1d, interp2d, RectBivariateSpline

plt.rcParams['figure.figsize'] = (6,5)
plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['legend.borderpad'] = 0.1
plt.rcParams['legend.labelspacing'] = 0.1
plt.rcParams['legend.handletextpad'] = 0.1
plt.rcParams['font.family'] = 'stixgeneral'
plt.rcParams['font.size'] = 16

temp = [1800,1900,2000,2100]

wave, ocl18w = np.loadtxt('ocl_trans_plots/1800_west_opticaldepth_level.txt', unpack = True)
wave, ocl18e = np.loadtxt('ocl_trans_plots/1800_east_opticaldepth_level.txt', unpack = True)
wave, ocl18p = np.loadtxt('ocl_trans_plots/1800_pole_opticaldepth_level.txt', unpack = True)
wave, ocl19w = np.loadtxt('ocl_trans_plots/1900_west_opticaldepth_level.txt', unpack = True)
wave, ocl19e = np.loadtxt('ocl_trans_plots/1900_east_opticaldepth_level.txt', unpack = True)
wave, ocl19p = np.loadtxt('ocl_trans_plots/1900_pole_opticaldepth_level.txt', unpack = True)
wave, ocl20w = np.loadtxt('ocl_trans_plots/2000_west_opticaldepth_level.txt', unpack = True)
wave, ocl20e = np.loadtxt('ocl_trans_plots/2000_east_opticaldepth_level.txt', unpack = True)
wave, ocl20p = np.loadtxt('ocl_trans_plots/2000_pole_opticaldepth_level.txt', unpack = True)
wave, ocl21w = np.loadtxt('ocl_trans_plots/2100_west_opticaldepth_level.txt', unpack = True)
wave, ocl21e = np.loadtxt('ocl_trans_plots/2100_east_opticaldepth_level.txt', unpack = True)
wave, ocl21p = np.loadtxt('ocl_trans_plots/2100_pole_opticaldepth_level.txt', unpack = True)

west = [ocl18w[3],ocl19w[3],ocl20w[3],ocl21w[3]]
east = [ocl18e[3],ocl19e[3],ocl20e[3],ocl21e[3]]
pole = [ocl18p[3],ocl19p[3],ocl20p[3],ocl21p[3]]


wave, ocl18w_min = np.loadtxt('ocl_trans_plots/1800_west_min_opticaldepth_level.txt', unpack = True)
wave, ocl18e_min = np.loadtxt('ocl_trans_plots/1800_east_min_opticaldepth_level.txt', unpack = True)
wave, ocl18p_min = np.loadtxt('ocl_trans_plots/1800_pole_min_opticaldepth_level.txt', unpack = True)
wave, ocl19w_min = np.loadtxt('ocl_trans_plots/1900_west_min_opticaldepth_level.txt', unpack = True)
wave, ocl19e_min = np.loadtxt('ocl_trans_plots/1900_east_min_opticaldepth_level.txt', unpack = True)
wave, ocl19p_min = np.loadtxt('ocl_trans_plots/1900_pole_min_opticaldepth_level.txt', unpack = True)
wave, ocl20w_min = np.loadtxt('ocl_trans_plots/2000_west_min_opticaldepth_level.txt', unpack = True)
wave, ocl20e_min = np.loadtxt('ocl_trans_plots/2000_east_min_opticaldepth_level.txt', unpack = True)
wave, ocl20p_min = np.loadtxt('ocl_trans_plots/2000_pole_min_opticaldepth_level.txt', unpack = True)
wave, ocl21w_min = np.loadtxt('ocl_trans_plots/2100_west_min_opticaldepth_level.txt', unpack = True)
wave, ocl21e_min = np.loadtxt('ocl_trans_plots/2100_east_min_opticaldepth_level.txt', unpack = True)
wave, ocl21p_min = np.loadtxt('ocl_trans_plots/2100_pole_min_opticaldepth_level.txt', unpack = True)

west_min = [ocl18w_min[3],ocl19w_min[3],ocl20w_min[3],ocl21w_min[3]]
east_min = [ocl18e_min[3],ocl19e_min[3],ocl20e_min[3],ocl21e_min[3]]
pole_min = [ocl18p_min[3],ocl19p_min[3],ocl20p_min[3],ocl21p_min[3]]

f,ax = plt.subplots(1,1)

ax.scatter(wave,ocl18w)
ax.plot(wave,ocl18w, c='r')
ax.scatter(wave,ocl18e)
ax.plot(wave,ocl18e, c='g')
ax.scatter(wave,ocl18p)
ax.plot(wave,ocl18p, c='b')
ax.scatter(wave,ocl18w_min)
ax.plot(wave,ocl18w_min, ls='--', c='r')
ax.scatter(wave,ocl18e_min)
ax.plot(wave,ocl18e_min, ls='--', c='g')
ax.scatter(wave,ocl18p_min)
ax.plot(wave,ocl18p, ls='--', c='b')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(min(wave),max(wave))
ax.set_ylim(1e-2,1e-7)
ax.set_ylabel(r'Pressure [bar]')
ax.set_xlabel(r'Wavelength [$\mu$m]')
plt.savefig('ocl_trans_plots/ocl_trend_1800.jpg', format = 'jpg', bbox_inches='tight')
plt.close()	

plt.rcParams['figure.figsize'] = (6,5)
plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['legend.borderpad'] = 0.1
plt.rcParams['legend.labelspacing'] = 0.1
plt.rcParams['legend.handletextpad'] = 0.1
plt.rcParams['font.family'] = 'stixgeneral'
plt.rcParams['font.size'] = 16
f,ax = plt.subplots(1,1)	

print wave[3]

ax.scatter(temp,west, c='#7570b3')
ax.plot(temp,west, c='#7570b3', lw=2)
ax.scatter(temp,east, c='#d95f02')
ax.plot(temp,east, c='#d95f02', lw=2)
ax.scatter(temp,pole, c='#1b9e77')
ax.plot(temp,pole, c='#1b9e77', lw=2)
#ax.scatter(temp,west_min, c='#7570b3')
#ax.plot(temp,west_min, c='#7570b3', ls='--')
#ax.scatter(temp,east_min, c='#d95f02')
#ax.plot(temp,east_min, c='#d95f02', ls='--')
#ax.scatter(temp,pole_min, c='#1b9e77')
#ax.plot(temp,pole_min, c='#1b9e77', ls='--')
ax.set_yscale('log')
ax.set_ylim(1e-2,1e-5)
ax.set_xlim(1750,2150)
ax.set_xticks((1800,1900,2000, 2100))
ax.set_xticklabels((1800,1900,2000, 2100))
ax.set_ylabel(r'Pressure [bar]')
ax.set_xlabel(r'Equilibrium Temperature [K]')
plt.savefig('ocl_trans_plots/ocl_trend.pdf', format = 'pdf', bbox_inches='tight')
plt.close()	


