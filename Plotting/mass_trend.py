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

temp, mass_w_max = np.loadtxt('mass_west_max.txt', unpack = True)
temp, mass_p_max = np.loadtxt('mass_pole_max.txt', unpack = True)
temp, mass_e_max = np.loadtxt('mass_east_max.txt', unpack = True)

#temp, mass_w_min = np.loadtxt('mass_west_min.txt', unpack = True)
#temp, mass_p_min = np.loadtxt('mass_pole_min.txt', unpack = True)
#temp, mass_e_min = np.loadtxt('mass_east_min.txt', unpack = True)

f, ax = plt.subplots(1,1)
ax.scatter(temp,mass_w_max,s=35, c='#7570b3')
ax.plot(temp,mass_w_max, lw=2, c='#7570b3', alpha = .75, label ='west limb, max')
ax.scatter(temp,mass_p_max,s=35, c='#1b9e77')
ax.plot(temp,mass_p_max, lw=2, c='#1b9e77', alpha = .75, label = 'pole, max')
ax.scatter(temp,mass_e_max,s=35, c='#d95f02')
ax.plot(temp,mass_e_max, lw=2, c='#d95f02', alpha = .75, label = 'east limb, max')
#ax.scatter(temp,mass_w_min,s=35, c='#7570b3')
#ax.plot(temp,mass_w_min, lw=2, c='#7570b3', alpha = .75, label ='west limb, min', ls='--')
#ax.scatter(temp,mass_p_min,s=35, c='#1b9e77')
#ax.plot(temp,mass_p_min, lw=2, c='#1b9e77', alpha = .75, label = 'pole, min', ls='--')
#ax.scatter(temp,mass_e_min,s=35, c='#d95f02')
#ax.plot(temp,mass_e_min, lw=2, c='#d95f02', alpha = .75, label = 'east limb, min', ls='--')
ax.set_yscale('log')
ax.set_ylim(1e-4,1e0)
ax.set_xlim(1750,2150)
ax.set_xticks((1800,1900,2000, 2100))
ax.set_xticklabels((1800,1900,2000, 2100))
#plt.legend(loc='best')
ax.set_xlabel(r'Equilibrium Temperature [K]')
ax.set_ylabel(r'Column Density [g cm$^{-2}$]')
plt.savefig('mass_trend_max.pdf', format = 'pdf', bbox_inches='tight')
plt.savefig('mass_trend_max.jpg', format = 'jpg', bbox_inches='tight')
plt.close()
