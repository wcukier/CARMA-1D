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

z, p_11, t_11, x, y = np.genfromtxt('../run/carma/t1100g100nc_m0.0.atm', skip_header = 1, skip_footer = 93, unpack = True)
z, p_12, t_12, x, y = np.genfromtxt('../run/carma/t1200g100nc_m0.0.atm', skip_header = 1, skip_footer = 93, unpack = True)
z, p_13, t_13, x, y = np.genfromtxt('../run/carma/t1300g100nc_m0.0.atm', skip_header = 1, skip_footer = 93, unpack = True)
z, p_14, t_14, x, y = np.genfromtxt('../run/carma/t1400g100nc_m0.0.atm', skip_header = 1, skip_footer = 93, unpack = True)
z, p_15, t_15, x, y = np.genfromtxt('../run/carma/t1500g100nc_m0.0.atm', skip_header = 1, skip_footer = 93, unpack = True)
z, p_16, t_16, x, y = np.genfromtxt('../run/carma/t1600g100nc_m0.0.atm', skip_header = 1, skip_footer = 93, unpack = True)


plt.plot(t_11,p_11*1e-6, c='k')
plt.plot(t_12,p_12*1e-6, c='k')
plt.plot(t_13,p_13*1e-6, c='k')
plt.plot(t_14,p_14*1e-6, c='k')
plt.plot(t_15,p_15*1e-6, c='k')
plt.plot(t_16,p_16*1e-6, c='k')

plt.plot((-32488)/(np.log10(5.936e-5/2)+np.log10(p_11*1e-6)+0.2*np.log10(p_11*1e-6)-14.88),p_11*1e-6, c='g', ls ='--', lw=2)
plt.plot(-74734.7/(np.log(p_11*1e-6+10**-7.08)-35.8027),p_11*1e-6, c='b', ls ='--', lw=2)
plt.plot((10**4)/(7.447-(0.42*np.log10(p_11*1e-6))),p_11*1e-6, c='y', ls ='--', lw=2)
plt.plot((-11382)/(np.log10(2.2e-7)+np.log10(p_11*1e-6)-7.611),p_11*1e-6, c='pink', ls ='--', lw=2)
plt.plot((-15873)/(np.log10(7.65e-8)+np.log10(p_11*1e-6)-12.812),p_11*1e-6, c='indigo', ls ='--', lw=2)
plt.plot((-45892.6)/(np.log10(4.937e-6/2)+np.log10(p_11*1e-6)-17.7),p_11*1e-6, c='turquoise', ls ='--', lw=2)
plt.plot((-20995)/(np.log10(5.78e-5)+np.log10(p_11*1e-6)-7.23),p_11*1e-6, c='orange', ls ='--', lw=2)
plt.plot((-20592)/(np.log10(8.87e-7)+np.log10(p_11*1e-6)-7.49),p_11*1e-6, c='r', ls ='--', lw=2)
plt.plot((-13889)/(np.log10(3.34e-6/2)+np.log10(p_11*1e-6)-8.55),p_11*1e-6, c='fuchsia', ls ='--', lw=2)
plt.yscale('log')
plt.ylim(max(p_11*1e-6), min(p_11*1e-6))
plt.xlabel(r"Temperature [K]", fontsize = 20)
plt.ylabel(r"Pressure [bar]", fontsize = 20)
plt.savefig('PT_profs.png', format = 'png', bbox_inches='tight')
plt.close()







