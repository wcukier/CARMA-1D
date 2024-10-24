#! /usr/bin/env python3.8
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


# Read in inputs from command line
opts, args = getopt.getopt(sys.argv[1:], "h")
file = args[0]
path = args[1]

# Check if folder where the plots are generated exists; if 
# not, make one. 
if not os.path.isdir(path):
	os.makedirs(path)
	
# Read in given output file and plot results. 
print('Read in output values file')
infile=open(file,'r')
line = infile.readline().split()

longitude = np.loadtxt('../run/carma/longitudes.txt')
ilong = len(longitude)

nz,ngroup,nelem,nbin,ngas,nstep,iskip = map(int,line)
r,ms,dr,rl,ru = zeros((nbin,ngroup)),zeros((nbin,ngroup)),zeros((nbin,ngroup)),zeros((nbin,ngroup)),zeros((nbin,ngroup))
p,t,z = zeros(nz),zeros(nz),zeros(nz)

ntime = int(nstep/iskip)

time = zeros(ntime)

M_KCL = 1.98
M_ZNS = 4.1
M_NA2S = 1.856
M_MNS = 3.3
M_CR = 7.19
M_FE = 7.874
M_MG2SIO4 = 3.21
M_TIO2 = 4.23
M_AL2O3 = 3.95

count = 0


for i in range(ngroup):
	for j in range(nbin):
		line = infile.readline().split()
		r[j,i] = float(line[2])
		ms[j,i] = float(line[3])
		dr[j,i] = float(line[4])
		rl[j,i] = float(line[5])
		ru[j,i] = float(line[6])
for i in range(nz):
	line = infile.readline().split()
	z[i] = float(line[1])
	p[i] = float(line[3])
	t[i] = float(line[4])

 
#PLOT SATURATION TEMPERATURE INSTEAD (as a function of pressure)
#WHERE THE CONDENSATION CURVE INTERSECTS THE T/P PROFILES (SHOW IF IT DOES THIS TWICE)
#need this P/T plots anyways


	
mg_temp = 1500.

temps_10 = np.loadtxt('../run/carma/equator_temperatures_1000.txt') #nz,nlong
mg_cond_temp_10 = np.zeros((ilong))

temps_10_top = []
p_top = []

for i in range(nz):
	if p[i]/1.e6 < 100.:
		p_top.append(p[i]/1.e6)
		temps_10_top.append(temps_10[i,:])	

p_top = np.array(p_top)
temps_10_top = np.array(temps_10_top)

for il, long in enumerate(longitude):
	f_temp = interp1d(temps_10[:,il],p/1.e6)
	mg_cond_temp_10[il] = f_temp(mg_temp)
	
temps_11 = np.loadtxt('../run/carma/equator_temperatures_1100.txt') #nz,nlong
mg_cond_temp_11 = np.zeros((ilong))

temps_11_top = []
p_top = []

for i in range(nz):
	if p[i]/1.e6 < 30.:
		p_top.append(p[i]/1.e6)
		temps_11_top.append(temps_10[i,:])	

p_top = np.array(p_top)
temps_11_top = np.array(temps_11_top)

print(temps_11_top[:,0],temps_11_top[:,2],temps_11_top[:,3])
mg_temp = 1300.
for il, long in enumerate(longitude):
	f_temp = interp1d(temps_11_top[:,il],p_top)
	mg_cond_temp_11[il] = f_temp(mg_temp)
	
mg_temp = 1500.
temps_12 = np.loadtxt('../run/carma/equator_temperatures_1200.txt') #nz,nlong
mg_cond_temp_12 = np.zeros((ilong))

for il, long in enumerate(longitude):
	f_temp = interp1d(temps_12[:,il],p/1.e6)
	mg_cond_temp_12[il] = f_temp(mg_temp)
	if mg_cond_temp_12[il] > 1.:
		mg_cond_temp_12[il] = mg_cond_temp_12[il-1]

temps_13 = np.loadtxt('../run/carma/equator_temperatures_1300.txt') #nz,nlong
mg_cond_temp_13 = np.zeros((ilong))

for il, long in enumerate(longitude):
	f_temp = interp1d(temps_13[:,il],p/1.e6)
	mg_cond_temp_13[il] = f_temp(mg_temp)

temps_14 = np.loadtxt('../run/carma/equator_temperatures_1400.txt') #nz,nlong
mg_cond_temp_14 = np.zeros((ilong))

for il, long in enumerate(longitude):
	f_temp = interp1d(temps_14[:,il],p/1.e6)
	mg_cond_temp_14[il] = f_temp(mg_temp)

temps_15 = np.loadtxt('../run/carma/equator_temperatures_1500.txt') #nz,nlong
mg_cond_temp_15 = np.zeros((ilong))

for il, long in enumerate(longitude):
	f_temp = interp1d(temps_15[:,il],p/1.e6)
	mg_cond_temp_15[il] = f_temp(mg_temp)

temps_16 = np.loadtxt('../run/carma/equator_temperatures_1600.txt') #nz,nlong
mg_cond_temp_16 = np.zeros((ilong))

for il, long in enumerate(longitude):
	f_temp = interp1d(temps_16[:,il],p/1.e6)
	mg_cond_temp_16[il] = f_temp(mg_temp)

temps_17 = np.loadtxt('../run/carma/equator_temperatures_1700.txt') #nz,nlong
mg_cond_temp_17 = np.zeros((ilong))

for il, long in enumerate(longitude):
	f_temp = interp1d(temps_17[:,il],p/1.e6)
	mg_cond_temp_17[il] = f_temp(mg_temp)
	
temps_18 = np.loadtxt('../run/carma/equator_temperatures_1800.txt') #nz,nlong
mg_cond_temp_18 = np.zeros((ilong))

for il, long in enumerate(longitude):
	f_temp = interp1d(temps_18[:,il],p/1.e6)
	mg_cond_temp_18[il] = f_temp(mg_temp)
 
temps_19 = np.loadtxt('../run/carma/equator_temperatures_1900.txt') #nz,nlong
mg_cond_temp_19 = np.zeros((ilong))

for il, long in enumerate(longitude):
	f_temp = interp1d(temps_19[:,il],p/1.e6)
	mg_cond_temp_19[il] = f_temp(mg_temp)
	
temps_20 = np.loadtxt('../run/carma/equator_temperatures_2000.txt') #nz,nlong
mg_cond_temp_20 = np.zeros((ilong))

for il, long in enumerate(longitude):
	f_temp = interp1d(temps_20[:,il],p/1.e6)
	if min(temps_20[:,il]) > mg_temp:
		mg_cond_temp_20[il] = min(p)/1.e6
	else:
		mg_cond_temp_20[il] = f_temp(mg_temp)
	
temps_21 = np.loadtxt('../run/carma/equator_temperatures_2100.txt') #nz,nlong
mg_cond_temp_21 = np.zeros((ilong))

for il, long in enumerate(longitude):
	f_temp = interp1d(temps_21[:,il],p/1.e6)
	if min(temps_21[:,il]) > mg_temp:
		mg_cond_temp_21[il] = min(p)/1.e6
	else:
		mg_cond_temp_21[il] = f_temp(mg_temp)




al_num_10 = np.loadtxt('1000_plots/num_density1000_al.txt')
al_num_10 = al_num_10.T #nz,nlong
cr_num_10 = np.loadtxt('1000_plots/num_density1000_cr.txt')
cr_num_10 = cr_num_10.T #nz,nlong
fe_num_10 = np.loadtxt('1000_plots/num_density1000_fe.txt')
fe_num_10 = fe_num_10.T #nz,nlong
kc_num_10 = np.loadtxt('1000_plots/num_density1000_kc.txt')
kc_num_10 = kc_num_10.T #nz,nlong
mg_num_10 = np.loadtxt('1000_plots/num_density1000_mg.txt')
mg_num_10 = mg_num_10.T #nz,nlong
mn_num_10 = np.loadtxt('1000_plots/num_density1000_mn.txt')
mn_num_10 = mn_num_10.T #nz,nlong
na_num_10 = np.loadtxt('1000_plots/num_density1000_na.txt')
na_num_10 = na_num_10.T #nz,nlong
ti_num_10 = np.loadtxt('1000_plots/num_density1000_ti.txt')
ti_num_10 = ti_num_10.T #nz,nlong
zn_num_10 = np.loadtxt('1000_plots/num_density1000_zn.txt')
zn_num_10 = zn_num_10.T #nz,nlong

al_mass_10 = np.loadtxt('1000_plots/mass_density1000_al.txt')
al_mass_10 = al_mass_10.T #nz,nlong
cr_mass_10 = np.loadtxt('1000_plots/mass_density1000_cr.txt')
cr_mass_10 = cr_mass_10.T #nz,nlong
fe_mass_10 = np.loadtxt('1000_plots/mass_density1000_fe.txt')
fe_mass_10 = fe_mass_10.T #nz,nlong
kc_mass_10 = np.loadtxt('1000_plots/mass_density1000_kc.txt')
kc_mass_10 = kc_mass_10.T #nz,nlong
mg_mass_10 = np.loadtxt('1000_plots/mass_density1000_mg.txt')
mg_mass_10 = mg_mass_10.T #nz,nlong
mn_mass_10 = np.loadtxt('1000_plots/mass_density1000_mn.txt')
mn_mass_10 = mn_mass_10.T #nz,nlong
na_mass_10 = np.loadtxt('1000_plots/mass_density1000_na.txt')
na_mass_10 = na_mass_10.T #nz,nlong
ti_mass_10 = np.loadtxt('1000_plots/mass_density1000_ti.txt')
ti_mass_10 = ti_mass_10.T #nz,nlong
zn_mass_10 = np.loadtxt('1000_plots/mass_density1000_zn.txt')
zn_mass_10 = zn_mass_10.T #nz,nlong

al_num_11 = np.loadtxt('1100_plots/num_density1100_al.txt')
al_num_11 = al_num_11.T #nz,nlong
cr_num_11 = np.loadtxt('1100_plots/num_density1100_cr.txt')
cr_num_11 = cr_num_11.T #nz,nlong
fe_num_11 = np.loadtxt('1100_plots/num_density1100_fe.txt')
fe_num_11 = fe_num_11.T #nz,nlong
kc_num_11 = np.loadtxt('1100_plots/num_density1100_kc.txt')
kc_num_11 = kc_num_11.T #nz,nlong
mg_num_11 = np.loadtxt('1100_plots/num_density1100_mg.txt')
mg_num_11 = mg_num_11.T #nz,nlong
mn_num_11 = np.loadtxt('1100_plots/num_density1100_mn.txt')
mn_num_11 = mn_num_11.T #nz,nlong
na_num_11 = np.loadtxt('1100_plots/num_density1100_na.txt')
na_num_11 = na_num_11.T #nz,nlong
ti_num_11 = np.loadtxt('1100_plots/num_density1100_ti.txt')
ti_num_11 = ti_num_11.T #nz,nlong
zn_num_11 = np.loadtxt('1100_plots/num_density1100_zn.txt')
zn_num_11 = zn_num_11.T #nz,nlong

al_mass_11 = np.loadtxt('1100_plots/mass_density1100_al.txt')
al_mass_11 = al_mass_11.T #nz,nlong
cr_mass_11 = np.loadtxt('1100_plots/mass_density1100_cr.txt')
cr_mass_11 = cr_mass_11.T #nz,nlong
fe_mass_11 = np.loadtxt('1100_plots/mass_density1100_fe.txt')
fe_mass_11 = fe_mass_11.T #nz,nlong
kc_mass_11 = np.loadtxt('1100_plots/mass_density1100_kc.txt')
kc_mass_11 = kc_mass_11.T #nz,nlong
mg_mass_11 = np.loadtxt('1100_plots/mass_density1100_mg.txt')
mg_mass_11 = mg_mass_11.T #nz,nlong
mn_mass_11 = np.loadtxt('1100_plots/mass_density1100_mn.txt')
mn_mass_11 = mn_mass_11.T #nz,nlong
na_mass_11 = np.loadtxt('1100_plots/mass_density1100_na.txt')
na_mass_11 = na_mass_11.T #nz,nlong
ti_mass_11 = np.loadtxt('1100_plots/mass_density1100_ti.txt')
ti_mass_11 = ti_mass_11.T #nz,nlong
zn_mass_11 = np.loadtxt('1100_plots/mass_density1100_zn.txt')
zn_mass_11 = zn_mass_11.T #nz,nlong

al_num_12 = np.loadtxt('1200_plots/num_density1200_al.txt')
al_num_12 = al_num_12.T #nz,nlong
cr_num_12 = np.loadtxt('1200_plots/num_density1200_cr.txt')
cr_num_12 = cr_num_12.T #nz,nlong
fe_num_12 = np.loadtxt('1200_plots/num_density1200_fe.txt')
fe_num_12 = fe_num_12.T #nz,nlong
kc_num_12 = np.loadtxt('1200_plots/num_density1200_kc.txt')
kc_num_12 = kc_num_12.T #nz,nlong
mg_num_12 = np.loadtxt('1200_plots/num_density1200_mg.txt')
mg_num_12 = mg_num_12.T #nz,nlong
mn_num_12 = np.loadtxt('1200_plots/num_density1200_mn.txt')
mn_num_12 = mn_num_12.T #nz,nlong
na_num_12 = np.loadtxt('1200_plots/num_density1200_na.txt')
na_num_12 = na_num_12.T #nz,nlong
ti_num_12 = np.loadtxt('1200_plots/num_density1200_ti.txt')
ti_num_12 = ti_num_12.T #nz,nlong
zn_num_12 = np.loadtxt('1200_plots/num_density1200_zn.txt')
zn_num_12 = zn_num_12.T #nz,nlong

al_mass_12 = np.loadtxt('1200_plots/mass_density1200_al.txt')
al_mass_12 = al_mass_12.T #nz,nlong
cr_mass_12 = np.loadtxt('1200_plots/mass_density1200_cr.txt')
cr_mass_12 = cr_mass_12.T #nz,nlong
fe_mass_12 = np.loadtxt('1200_plots/mass_density1200_fe.txt')
fe_mass_12 = fe_mass_12.T #nz,nlong
kc_mass_12 = np.loadtxt('1200_plots/mass_density1200_kc.txt')
kc_mass_12 = kc_mass_12.T #nz,nlong
mg_mass_12 = np.loadtxt('1200_plots/mass_density1200_mg.txt')
mg_mass_12 = mg_mass_12.T #nz,nlong
mn_mass_12 = np.loadtxt('1200_plots/mass_density1200_mn.txt')
mn_mass_12 = mn_mass_12.T #nz,nlong
na_mass_12 = np.loadtxt('1200_plots/mass_density1200_na.txt')
na_mass_12 = na_mass_12.T #nz,nlong
ti_mass_12 = np.loadtxt('1200_plots/mass_density1200_ti.txt')
ti_mass_12 = ti_mass_12.T #nz,nlong
zn_mass_12 = np.loadtxt('1200_plots/mass_density1200_zn.txt')
zn_mass_12 = zn_mass_12.T #nz,nlong

al_num_13 = np.loadtxt('1300_plots/num_density1300_al.txt')
al_num_13 = al_num_13.T #nz,nlong
cr_num_13 = np.loadtxt('1300_plots/num_density1300_cr.txt')
cr_num_13 = cr_num_13.T #nz,nlong
fe_num_13 = np.loadtxt('1300_plots/num_density1300_fe.txt')
fe_num_13 = fe_num_13.T #nz,nlong
kc_num_13 = np.loadtxt('1300_plots/num_density1300_kc.txt')
kc_num_13 = kc_num_13.T #nz,nlong
mg_num_13 = np.loadtxt('1300_plots/num_density1300_mg.txt')
mg_num_13 = mg_num_13.T #nz,nlong
mn_num_13 = np.loadtxt('1300_plots/num_density1300_mn.txt')
mn_num_13 = mn_num_13.T #nz,nlong
na_num_13 = np.loadtxt('1300_plots/num_density1300_na.txt')
na_num_13 = na_num_13.T #nz,nlong
ti_num_13 = np.loadtxt('1300_plots/num_density1300_ti.txt')
ti_num_13 = ti_num_13.T #nz,nlong
zn_num_13 = np.loadtxt('1300_plots/num_density1300_zn.txt')
zn_num_13 = zn_num_13.T #nz,nlong

al_mass_13 = np.loadtxt('1300_plots/mass_density1300_al.txt')
al_mass_13 = al_mass_13.T #nz,nlong
cr_mass_13 = np.loadtxt('1300_plots/mass_density1300_cr.txt')
cr_mass_13 = cr_mass_13.T #nz,nlong
fe_mass_13 = np.loadtxt('1300_plots/mass_density1300_fe.txt')
fe_mass_13 = fe_mass_13.T #nz,nlong
kc_mass_13 = np.loadtxt('1300_plots/mass_density1300_kc.txt')
kc_mass_13 = kc_mass_13.T #nz,nlong
mg_mass_13 = np.loadtxt('1300_plots/mass_density1300_mg.txt')
mg_mass_13 = mg_mass_13.T #nz,nlong
mn_mass_13 = np.loadtxt('1300_plots/mass_density1300_mn.txt')
mn_mass_13 = mn_mass_13.T #nz,nlong
na_mass_13 = np.loadtxt('1300_plots/mass_density1300_na.txt')
na_mass_13 = na_mass_13.T #nz,nlong
ti_mass_13 = np.loadtxt('1300_plots/mass_density1300_ti.txt')
ti_mass_13 = ti_mass_13.T #nz,nlong
zn_mass_13 = np.loadtxt('1300_plots/mass_density1300_zn.txt')
zn_mass_13 = zn_mass_13.T #nz,nlong

al_num_14 = np.loadtxt('1400_plots/num_density1400_al.txt')
al_num_14 = al_num_14.T #nz,nlong
cr_num_14 = np.loadtxt('1400_plots/num_density1400_cr.txt')
cr_num_14 = cr_num_14.T #nz,nlong
fe_num_14 = np.loadtxt('1400_plots/num_density1400_fe.txt')
fe_num_14 = fe_num_14.T #nz,nlong
kc_num_14 = np.loadtxt('1400_plots/num_density1400_kc.txt')
kc_num_14 = kc_num_14.T #nz,nlong
mg_num_14 = np.loadtxt('1400_plots/num_density1400_mg.txt')
mg_num_14 = mg_num_14.T #nz,nlong
mn_num_14 = np.loadtxt('1400_plots/num_density1400_mn.txt')
mn_num_14 = mn_num_14.T #nz,nlong
na_num_14 = np.loadtxt('1400_plots/num_density1400_na.txt')
na_num_14 = na_num_14.T #nz,nlong
ti_num_14 = np.loadtxt('1400_plots/num_density1400_ti.txt')
ti_num_14 = ti_num_14.T #nz,nlong
zn_num_14 = np.loadtxt('1400_plots/num_density1400_zn.txt')
zn_num_14 = zn_num_14.T #nz,nlong

al_mass_14 = np.loadtxt('1400_plots/mass_density1400_al.txt')
al_mass_14 = al_mass_14.T #nz,nlong
cr_mass_14 = np.loadtxt('1400_plots/mass_density1400_cr.txt')
cr_mass_14 = cr_mass_14.T #nz,nlong
fe_mass_14 = np.loadtxt('1400_plots/mass_density1400_fe.txt')
fe_mass_14 = fe_mass_14.T #nz,nlong
kc_mass_14 = np.loadtxt('1400_plots/mass_density1400_kc.txt')
kc_mass_14 = kc_mass_14.T #nz,nlong
mg_mass_14 = np.loadtxt('1400_plots/mass_density1400_mg.txt')
mg_mass_14 = mg_mass_14.T #nz,nlong
mn_mass_14 = np.loadtxt('1400_plots/mass_density1400_mn.txt')
mn_mass_14 = mn_mass_14.T #nz,nlong
na_mass_14 = np.loadtxt('1400_plots/mass_density1400_na.txt')
na_mass_14 = na_mass_14.T #nz,nlong
ti_mass_14 = np.loadtxt('1400_plots/mass_density1400_ti.txt')
ti_mass_14 = ti_mass_14.T #nz,nlong
zn_mass_14 = np.loadtxt('1400_plots/mass_density1400_zn.txt')
zn_mass_14 = zn_mass_14.T #nz,nlong

al_num_15 = np.loadtxt('1500_plots/num_density1500_al.txt')
al_num_15 = al_num_15.T #nz,nlong
cr_num_15 = np.loadtxt('1500_plots/num_density1500_cr.txt')
cr_num_15 = cr_num_15.T #nz,nlong
fe_num_15 = np.loadtxt('1500_plots/num_density1500_fe.txt')
fe_num_15 = fe_num_15.T #nz,nlong
kc_num_15 = np.loadtxt('1500_plots/num_density1500_kc.txt')
kc_num_15 = kc_num_15.T #nz,nlong
mg_num_15 = np.loadtxt('1500_plots/num_density1500_mg.txt')
mg_num_15 = mg_num_15.T #nz,nlong
mn_num_15 = np.loadtxt('1500_plots/num_density1500_mn.txt')
mn_num_15 = mn_num_15.T #nz,nlong
na_num_15 = np.loadtxt('1500_plots/num_density1500_na.txt')
na_num_15 = na_num_15.T #nz,nlong
ti_num_15 = np.loadtxt('1500_plots/num_density1500_ti.txt')
ti_num_15 = ti_num_15.T #nz,nlong
zn_num_15 = np.loadtxt('1500_plots/num_density1500_zn.txt')
zn_num_15 = zn_num_15.T #nz,nlong

al_mass_15 = np.loadtxt('1500_plots/mass_density1500_al.txt')
al_mass_15 = al_mass_15.T #nz,nlong
cr_mass_15 = np.loadtxt('1500_plots/mass_density1500_cr.txt')
cr_mass_15 = cr_mass_15.T #nz,nlong
fe_mass_15 = np.loadtxt('1500_plots/mass_density1500_fe.txt')
fe_mass_15 = fe_mass_15.T #nz,nlong
kc_mass_15 = np.loadtxt('1500_plots/mass_density1500_kc.txt')
kc_mass_15 = kc_mass_15.T #nz,nlong
mg_mass_15 = np.loadtxt('1500_plots/mass_density1500_mg.txt')
mg_mass_15 = mg_mass_15.T #nz,nlong
mn_mass_15 = np.loadtxt('1500_plots/mass_density1500_mn.txt')
mn_mass_15 = mn_mass_15.T #nz,nlong
na_mass_15 = np.loadtxt('1500_plots/mass_density1500_na.txt')
na_mass_15 = na_mass_15.T #nz,nlong
ti_mass_15 = np.loadtxt('1500_plots/mass_density1500_ti.txt')
ti_mass_15 = ti_mass_15.T #nz,nlong
zn_mass_15 = np.loadtxt('1500_plots/mass_density1500_zn.txt')
zn_mass_15 = zn_mass_15.T #nz,nlong

al_num_16 = np.loadtxt('1600_plots/num_density1600_al.txt')
al_num_16 = al_num_16.T #nz,nlong
cr_num_16 = np.loadtxt('1600_plots/num_density1600_cr.txt')
cr_num_16 = cr_num_16.T #nz,nlong
fe_num_16 = np.loadtxt('1600_plots/num_density1600_fe.txt')
fe_num_16 = fe_num_16.T #nz,nlong
kc_num_16 = np.loadtxt('1600_plots/num_density1600_kc.txt')
kc_num_16 = kc_num_16.T #nz,nlong
mg_num_16 = np.loadtxt('1600_plots/num_density1600_mg.txt')
mg_num_16 = mg_num_16.T #nz,nlong
mn_num_16 = np.loadtxt('1600_plots/num_density1600_mn.txt')
mn_num_16 = mn_num_16.T #nz,nlong
na_num_16 = np.loadtxt('1600_plots/num_density1600_na.txt')
na_num_16 = na_num_16.T #nz,nlong
ti_num_16 = np.loadtxt('1600_plots/num_density1600_ti.txt')
ti_num_16 = ti_num_16.T #nz,nlong
zn_num_16 = np.loadtxt('1600_plots/num_density1600_zn.txt')
zn_num_16 = zn_num_16.T #nz,nlong

al_mass_16 = np.loadtxt('1600_plots/mass_density1600_al.txt')
al_mass_16 = al_mass_16.T #nz,nlong
cr_mass_16 = np.loadtxt('1600_plots/mass_density1600_cr.txt')
cr_mass_16 = cr_mass_16.T #nz,nlong
fe_mass_16 = np.loadtxt('1600_plots/mass_density1600_fe.txt')
fe_mass_16 = fe_mass_16.T #nz,nlong
kc_mass_16 = np.loadtxt('1600_plots/mass_density1600_kc.txt')
kc_mass_16 = kc_mass_16.T #nz,nlong
mg_mass_16 = np.loadtxt('1600_plots/mass_density1600_mg.txt')
mg_mass_16 = mg_mass_16.T #nz,nlong
mn_mass_16 = np.loadtxt('1600_plots/mass_density1600_mn.txt')
mn_mass_16 = mn_mass_16.T #nz,nlong
na_mass_16 = np.loadtxt('1600_plots/mass_density1600_na.txt')
na_mass_16 = na_mass_16.T #nz,nlong
ti_mass_16 = np.loadtxt('1600_plots/mass_density1600_ti.txt')
ti_mass_16 = ti_mass_16.T #nz,nlong
zn_mass_16 = np.loadtxt('1600_plots/mass_density1600_zn.txt')
zn_mass_16 = zn_mass_16.T #nz,nlong

al_num_17 = np.loadtxt('1700_plots/num_density1700_al.txt')
al_num_17 = al_num_17.T #nz,nlong
cr_num_17 = np.loadtxt('1700_plots/num_density1700_cr.txt')
cr_num_17 = cr_num_17.T #nz,nlong
fe_num_17 = np.loadtxt('1700_plots/num_density1700_fe.txt')
fe_num_17 = fe_num_17.T #nz,nlong
kc_num_17 = np.loadtxt('1700_plots/num_density1700_kc.txt')
kc_num_17 = kc_num_17.T #nz,nlong
mg_num_17 = np.loadtxt('1700_plots/num_density1700_mg.txt')
mg_num_17 = mg_num_17.T #nz,nlong
mn_num_17 = np.loadtxt('1700_plots/num_density1700_mn.txt')
mn_num_17 = mn_num_17.T #nz,nlong
na_num_17 = np.loadtxt('1700_plots/num_density1700_na.txt')
na_num_17 = na_num_17.T #nz,nlong
ti_num_17 = np.loadtxt('1700_plots/num_density1700_ti.txt')
ti_num_17 = ti_num_17.T #nz,nlong
zn_num_17 = np.loadtxt('1700_plots/num_density1700_zn.txt')
zn_num_17 = zn_num_17.T #nz,nlong

al_mass_17 = np.loadtxt('1700_plots/mass_density1700_al.txt')
al_mass_17 = al_mass_17.T #nz,nlong
cr_mass_17 = np.loadtxt('1700_plots/mass_density1700_cr.txt')
cr_mass_17 = cr_mass_17.T #nz,nlong
fe_mass_17 = np.loadtxt('1700_plots/mass_density1700_fe.txt')
fe_mass_17 = fe_mass_17.T #nz,nlong
kc_mass_17 = np.loadtxt('1700_plots/mass_density1700_kc.txt')
kc_mass_17 = kc_mass_17.T #nz,nlong
mg_mass_17 = np.loadtxt('1700_plots/mass_density1700_mg.txt')
mg_mass_17 = mg_mass_17.T #nz,nlong
mn_mass_17 = np.loadtxt('1700_plots/mass_density1700_mn.txt')
mn_mass_17 = mn_mass_17.T #nz,nlong
na_mass_17 = np.loadtxt('1700_plots/mass_density1700_na.txt')
na_mass_17 = na_mass_17.T #nz,nlong
ti_mass_17 = np.loadtxt('1700_plots/mass_density1700_ti.txt')
ti_mass_17 = ti_mass_17.T #nz,nlong
zn_mass_17 = np.loadtxt('1700_plots/mass_density1700_zn.txt')
zn_mass_17 = zn_mass_17.T #nz,nlong

al_num_18 = np.loadtxt('1800_plots/num_density1800_al.txt')
al_num_18 = al_num_18.T #nz,nlong
cr_num_18 = np.loadtxt('1800_plots/num_density1800_cr.txt')
cr_num_18 = cr_num_18.T #nz,nlong
fe_num_18 = np.loadtxt('1800_plots/num_density1800_fe.txt')
fe_num_18 = fe_num_18.T #nz,nlong
kc_num_18 = np.loadtxt('1800_plots/num_density1800_kc.txt')
kc_num_18 = kc_num_18.T #nz,nlong
mg_num_18 = np.loadtxt('1800_plots/num_density1800_mg.txt')
mg_num_18 = mg_num_18.T #nz,nlong
mn_num_18 = np.loadtxt('1800_plots/num_density1800_mn.txt')
mn_num_18 = mn_num_18.T #nz,nlong
na_num_18 = np.loadtxt('1800_plots/num_density1800_na.txt')
na_num_18 = na_num_18.T #nz,nlong
ti_num_18 = np.loadtxt('1800_plots/num_density1800_ti.txt')
ti_num_18 = ti_num_18.T #nz,nlong
zn_num_18 = np.loadtxt('1800_plots/num_density1800_zn.txt')
zn_num_18 = zn_num_18.T #nz,nlong

al_mass_18 = np.loadtxt('1800_plots/mass_density1800_al.txt')
al_mass_18 = al_mass_18.T #nz,nlong
cr_mass_18 = np.loadtxt('1800_plots/mass_density1800_cr.txt')
cr_mass_18 = cr_mass_18.T #nz,nlong
fe_mass_18 = np.loadtxt('1800_plots/mass_density1800_fe.txt')
fe_mass_18 = fe_mass_18.T #nz,nlong
kc_mass_18 = np.loadtxt('1800_plots/mass_density1800_kc.txt')
kc_mass_18 = kc_mass_18.T #nz,nlong
mg_mass_18 = np.loadtxt('1800_plots/mass_density1800_mg.txt')
mg_mass_18 = mg_mass_18.T #nz,nlong
mn_mass_18 = np.loadtxt('1800_plots/mass_density1800_mn.txt')
mn_mass_18 = mn_mass_18.T #nz,nlong
na_mass_18 = np.loadtxt('1800_plots/mass_density1800_na.txt')
na_mass_18 = na_mass_18.T #nz,nlong
ti_mass_18 = np.loadtxt('1800_plots/mass_density1800_ti.txt')
ti_mass_18 = ti_mass_18.T #nz,nlong
zn_mass_18 = np.loadtxt('1800_plots/mass_density1800_zn.txt')
zn_mass_18 = zn_mass_18.T #nz,nlong

al_num_19 = np.loadtxt('1900_plots/num_density1900_al.txt')
al_num_19 = al_num_19.T #nz,nlong
cr_num_19 = np.loadtxt('1900_plots/num_density1900_cr.txt')
cr_num_19 = cr_num_19.T #nz,nlong
fe_num_19 = np.loadtxt('1900_plots/num_density1900_fe.txt')
fe_num_19 = fe_num_19.T #nz,nlong
kc_num_19 = np.loadtxt('1900_plots/num_density1900_kc.txt')
kc_num_19 = kc_num_19.T #nz,nlong
mg_num_19 = np.loadtxt('1900_plots/num_density1900_mg.txt')
mg_num_19 = mg_num_19.T #nz,nlong
mn_num_19 = np.loadtxt('1900_plots/num_density1900_mn.txt')
mn_num_19 = mn_num_19.T #nz,nlong
na_num_19 = np.loadtxt('1900_plots/num_density1900_na.txt')
na_num_19 = na_num_19.T #nz,nlong
ti_num_19 = np.loadtxt('1900_plots/num_density1900_ti.txt')
ti_num_19 = ti_num_19.T #nz,nlong
zn_num_19 = np.loadtxt('1900_plots/num_density1900_zn.txt')
zn_num_19 = zn_num_19.T #nz,nlong

al_mass_19 = np.loadtxt('1900_plots/mass_density1900_al.txt')
al_mass_19 = al_mass_19.T #nz,nlong
cr_mass_19 = np.loadtxt('1900_plots/mass_density1900_cr.txt')
cr_mass_19 = cr_mass_19.T #nz,nlong
fe_mass_19 = np.loadtxt('1900_plots/mass_density1900_fe.txt')
fe_mass_19 = fe_mass_19.T #nz,nlong
kc_mass_19 = np.loadtxt('1900_plots/mass_density1900_kc.txt')
kc_mass_19 = kc_mass_19.T #nz,nlong
mg_mass_19 = np.loadtxt('1900_plots/mass_density1900_mg.txt')
mg_mass_19 = mg_mass_19.T #nz,nlong
mn_mass_19 = np.loadtxt('1900_plots/mass_density1900_mn.txt')
mn_mass_19 = mn_mass_19.T #nz,nlong
na_mass_19 = np.loadtxt('1900_plots/mass_density1900_na.txt')
na_mass_19 = na_mass_19.T #nz,nlong
ti_mass_19 = np.loadtxt('1900_plots/mass_density1900_ti.txt')
ti_mass_19 = ti_mass_19.T #nz,nlong
zn_mass_19 = np.loadtxt('1900_plots/mass_density1900_zn.txt')
zn_mass_19 = zn_mass_19.T #nz,nlong

al_num_20 = np.loadtxt('2000_plots/num_density2000_al.txt')
al_num_20 = al_num_20.T #nz,nlong
cr_num_20 = np.loadtxt('2000_plots/num_density2000_cr.txt')
cr_num_20 = cr_num_20.T #nz,nlong
fe_num_20 = np.loadtxt('2000_plots/num_density2000_fe.txt')
fe_num_20 = fe_num_20.T #nz,nlong
kc_num_20 = np.loadtxt('2000_plots/num_density2000_kc.txt')
kc_num_20 = kc_num_20.T #nz,nlong
mg_num_20 = np.loadtxt('2000_plots/num_density2000_mg.txt')
mg_num_20 = mg_num_20.T #nz,nlong
mn_num_20 = np.loadtxt('2000_plots/num_density2000_mn.txt')
mn_num_20 = mn_num_20.T #nz,nlong
na_num_20 = np.loadtxt('2000_plots/num_density2000_na.txt')
na_num_20 = na_num_20.T #nz,nlong
ti_num_20 = np.loadtxt('2000_plots/num_density2000_ti.txt')
ti_num_20 = ti_num_20.T #nz,nlong
zn_num_20 = np.loadtxt('2000_plots/num_density2000_zn.txt')
zn_num_20 = zn_num_20.T #nz,nlong

al_mass_20 = np.loadtxt('2000_plots/mass_density2000_al.txt')
al_mass_20 = al_mass_20.T #nz,nlong
cr_mass_20 = np.loadtxt('2000_plots/mass_density2000_cr.txt')
cr_mass_20 = cr_mass_20.T #nz,nlong
fe_mass_20 = np.loadtxt('2000_plots/mass_density2000_fe.txt')
fe_mass_20 = fe_mass_20.T #nz,nlong
kc_mass_20 = np.loadtxt('2000_plots/mass_density2000_kc.txt')
kc_mass_20 = kc_mass_20.T #nz,nlong
mg_mass_20 = np.loadtxt('2000_plots/mass_density2000_mg.txt')
mg_mass_20 = mg_mass_20.T #nz,nlong
mn_mass_20 = np.loadtxt('2000_plots/mass_density2000_mn.txt')
mn_mass_20 = mn_mass_20.T #nz,nlong
na_mass_20 = np.loadtxt('2000_plots/mass_density2000_na.txt')
na_mass_20 = na_mass_20.T #nz,nlong
ti_mass_20 = np.loadtxt('2000_plots/mass_density2000_ti.txt')
ti_mass_20 = ti_mass_20.T #nz,nlong
zn_mass_20 = np.loadtxt('2000_plots/mass_density2000_zn.txt')
zn_mass_20 = zn_mass_20.T #nz,nlong

al_num_21 = np.loadtxt('2100_plots/num_density2100_al.txt')
al_num_21 = al_num_21.T #nz,nlong
cr_num_21 = np.loadtxt('2100_plots/num_density2100_cr.txt')
cr_num_21 = cr_num_21.T #nz,nlong
fe_num_21 = np.loadtxt('2100_plots/num_density2100_fe.txt')
fe_num_21 = fe_num_21.T #nz,nlong
kc_num_21 = np.loadtxt('2100_plots/num_density2100_kc.txt')
kc_num_21 = kc_num_21.T #nz,nlong
mg_num_21 = np.loadtxt('2100_plots/num_density2100_mg.txt')
mg_num_21 = mg_num_21.T #nz,nlong
mn_num_21 = np.loadtxt('2100_plots/num_density2100_mn.txt')
mn_num_21 = mn_num_21.T #nz,nlong
na_num_21 = np.loadtxt('2100_plots/num_density2100_na.txt')
na_num_21 = na_num_21.T #nz,nlong
ti_num_21 = np.loadtxt('2100_plots/num_density2100_ti.txt')
ti_num_21 = ti_num_21.T #nz,nlong
zn_num_21 = np.loadtxt('2100_plots/num_density2100_zn.txt')
zn_num_21 = zn_num_21.T #nz,nlong

al_mass_21 = np.loadtxt('2100_plots/mass_density2100_al.txt')
al_mass_21 = al_mass_21.T #nz,nlong
cr_mass_21 = np.loadtxt('2100_plots/mass_density2100_cr.txt')
cr_mass_21 = cr_mass_21.T #nz,nlong
fe_mass_21 = np.loadtxt('2100_plots/mass_density2100_fe.txt')
fe_mass_21 = fe_mass_21.T #nz,nlong
kc_mass_21 = np.loadtxt('2100_plots/mass_density2100_kc.txt')
kc_mass_21 = kc_mass_21.T #nz,nlong
mg_mass_21 = np.loadtxt('2100_plots/mass_density2100_mg.txt')
mg_mass_21 = mg_mass_21.T #nz,nlong
mn_mass_21 = np.loadtxt('2100_plots/mass_density2100_mn.txt')
mn_mass_21 = mn_mass_21.T #nz,nlong
na_mass_21 = np.loadtxt('2100_plots/mass_density2100_na.txt')
na_mass_21 = na_mass_21.T #nz,nlong
ti_mass_21 = np.loadtxt('2100_plots/mass_density2100_ti.txt')
ti_mass_21 = ti_mass_21.T #nz,nlong
zn_mass_21 = np.loadtxt('2100_plots/mass_density2100_zn.txt')
zn_mass_21 = zn_mass_21.T #nz,nlong


#########################################################################################################################################################################
# Read in given output file and plot results. 
print('PLOTTING!!!!')


levels_num = np.array([1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8])/1e7
levels_mass = np.array([1e-16,1e-14,1e-12,1e-10,1e-8,1e-6,1e-4,1e-2])/1e7

X1, Y1 = np.meshgrid(longitude,p/1.e6)
X2, Y2 = np.meshgrid(longitude,p/1.e6)
X3, Y3 = np.meshgrid(longitude,p/1.e6)
X4, Y4 = np.meshgrid(longitude,p/1.e6)
X5, Y5 = np.meshgrid(longitude,p/1.e6)
X6, Y6 = np.meshgrid(longitude,p/1.e6)
X7, Y7 = np.meshgrid(longitude,p/1.e6)
X8, Y8 = np.meshgrid(longitude,p/1.e6)
X9, Y9 = np.meshgrid(longitude,p/1.e6)
X10, Y10 = np.meshgrid(longitude,p/1.e6)
X11, Y11 = np.meshgrid(longitude,p/1.e6)

fig = plt.figure()
fig, ((ax1,ax2),(ax3,ax4),(ax5,ax6),(ax7,ax8),(ax9,ax10),(ax11,ax12)) = plt.subplots(nrows = 6, ncols = 2,sharey=True, sharex= True, figsize=(15,32))

fontsized = 24
plt.rcParams['font.size'] = fontsized
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'stixgeneral'

plt.rcParams['axes.linewidth'] = 1



im2 = ax1.contourf(X2,Y2,al_mass_10[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Greens', alpha=0.75)
im1 = ax1.contourf(X1,Y1,ti_mass_10[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Blues', alpha=0.75)
#im4 = ax1.contourf(X4,Y4,mg_mass_10[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
im3 = ax1.contourf(X3,Y3,fe_mass_10[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Reds', alpha=0.75)
im5 = ax1.contourf(X5,Y5,cr_mass_10[:,:],norm=LogNorm(),levels=levels_mass*10, cmap='Oranges', alpha=0.75)
im6 = ax1.contourf(X6,Y6,mn_mass_10[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Greys', alpha=0.75)
ax1.set_yscale('log')
ax1.set_ylim(100.,min(p)/1.e6)
ax1.set_xlim(min(longitude),max(longitude))
ax1.text(100,1e-5,'1000 K')


im2 = ax2.contourf(X2,Y2,al_mass_11[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Greens', alpha=0.75)
im1 = ax2.contourf(X1,Y1,ti_mass_11[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Blues', alpha=0.75)
#im4 = ax2.contourf(X4,Y4,mg_mass_11[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
im3 = ax2.contourf(X3,Y3,fe_mass_11[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Reds', alpha=0.75)
im5 = ax2.contourf(X5,Y5,cr_mass_11[:,:],norm=LogNorm(),levels=levels_mass*10, cmap='Oranges', alpha=0.75)
ax2.set_yscale('log')
ax2.set_ylim(100.,min(p)/1.e6)
ax2.set_xlim(min(longitude),max(longitude))
ax2.text(100,1e-5,'1100 K')


im2 = ax3.contourf(X2,Y2,al_mass_12[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Greens', alpha=0.75)
im1 = ax3.contourf(X1,Y1,ti_mass_12[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Blues', alpha=0.75)
#im4 = ax3.contourf(X4,Y4,mg_mass_12[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
im3 = ax3.contourf(X3,Y3,fe_mass_12[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Reds', alpha=0.75)
im5 = ax3.contourf(X5,Y5,cr_mass_12[:,:],norm=LogNorm(),levels=levels_mass*10, cmap='Oranges', alpha=0.75)
ax3.set_yscale('log')
ax3.set_ylim(100.,min(p)/1.e6)
ax3.set_xlim(min(longitude),max(longitude))
ax3.text(100,1e-5,'1200 K')


im2 = ax4.contourf(X2,Y2,al_mass_13[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Greens', alpha=0.75)
im1 = ax4.contourf(X1,Y1,ti_mass_13[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Blues', alpha=0.75)
#im4 = ax4.contourf(X4,Y4,mg_mass_13[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
im3 = ax4.contourf(X3,Y3,fe_mass_13[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Reds', alpha=0.75)
im5 = ax4.contourf(X5,Y5,cr_mass_13[:,:],norm=LogNorm(),levels=levels_mass*10, cmap='Oranges', alpha=0.75)
ax4.set_yscale('log')
ax4.set_ylim(100.,min(p)/1.e6)
ax4.set_xlim(min(longitude),max(longitude))
ax4.text(100,1e-5,'1300 K')

im2 = ax5.contourf(X2,Y2,al_mass_14[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Greens', alpha=0.75)
im1 = ax5.contourf(X1,Y1,ti_mass_14[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Blues', alpha=0.75)
#im4 = ax5.contourf(X4,Y4,mg_mass_14[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
im3 = ax5.contourf(X3,Y3,fe_mass_14[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Reds', alpha=0.75)
im5 = ax5.contourf(X5,Y5,cr_mass_14[:,:],norm=LogNorm(),levels=levels_mass*10, cmap='Oranges', alpha=0.75)
ax5.set_yscale('log')
ax5.set_ylim(100.,min(p)/1.e6)
ax5.set_xlim(min(longitude),max(longitude))
ax5.text(100,1e-5,'1400 K')

im2 = ax6.contourf(X2,Y2,al_mass_15[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Greens', alpha=0.75)
im1 = ax6.contourf(X1,Y1,ti_mass_15[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Blues', alpha=0.75)
#im4 = ax6.contourf(X4,Y4,mg_mass_15[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
im3 = ax6.contourf(X3,Y3,fe_mass_15[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Reds', alpha=0.75)
im5 = ax6.contourf(X5,Y5,cr_mass_15[:,:],norm=LogNorm(),levels=levels_mass*10, cmap='Oranges', alpha=0.75)
ax6.set_yscale('log')
ax6.set_ylim(100.,min(p)/1.e6)
ax6.set_xlim(min(longitude),max(longitude))
ax6.text(100,1e-5,'1500 K')

im2 = ax7.contourf(X2,Y2,al_mass_16[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Greens', alpha=0.75)
im1 = ax7.contourf(X1,Y1,ti_mass_16[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Blues', alpha=0.75)
#im4 = ax7.contourf(X4,Y4,mg_mass_16[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
im3 = ax7.contourf(X3,Y3,fe_mass_16[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Reds', alpha=0.75)
im5 = ax7.contourf(X5,Y5,cr_mass_16[:,:],norm=LogNorm(),levels=levels_mass*10, cmap='Oranges', alpha=0.75)
ax7.set_yscale('log')
ax7.set_ylim(100.,min(p)/1.e6)
ax7.set_xlim(min(longitude),max(longitude))
ax7.text(100,1e-5,'1600 K')

im2 = ax8.contourf(X2,Y2,al_mass_17[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Greens', alpha=0.75)
im1 = ax8.contourf(X1,Y1,ti_mass_17[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Blues', alpha=0.75)
#im4 = ax8.contourf(X4,Y4,mg_mass_17[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
im3 = ax8.contourf(X3,Y3,fe_mass_17[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Reds', alpha=0.75)
im5 = ax8.contourf(X5,Y5,cr_mass_17[:,:],norm=LogNorm(),levels=levels_mass*10, cmap='Oranges', alpha=0.75)
ax8.set_yscale('log')
ax8.set_ylim(100.,min(p)/1.e6)
ax8.set_xlim(min(longitude),max(longitude))
ax8.text(100,1e-5,'1700 K')

im2 = ax9.contourf(X2,Y2,al_mass_18[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Greens', alpha=0.75)
im1 = ax9.contourf(X1,Y1,ti_mass_18[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Blues', alpha=0.75)
#im4 = ax9.contourf(X4,Y4,mg_mass_18[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
im3 = ax9.contourf(X3,Y3,fe_mass_18[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Reds', alpha=0.75)
im5 = ax9.contourf(X5,Y5,cr_mass_18[:,:],norm=LogNorm(),levels=levels_mass*10, cmap='Oranges', alpha=0.75)
ax9.set_yscale('log')
ax9.set_ylim(100.,min(p)/1.e6)
ax9.set_xlim(min(longitude),max(longitude))
ax9.text(100,1e-5,'1800 K')

im2 = ax10.contourf(X2,Y2,al_mass_19[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Greens', alpha=0.75)
im1 = ax10.contourf(X1,Y1,ti_mass_19[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Blues', alpha=0.75)
#im4 = ax10.contourf(X4,Y4,mg_mass_19[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
im3 = ax10.contourf(X3,Y3,fe_mass_19[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Reds', alpha=0.75)
im5 = ax10.contourf(X5,Y5,cr_mass_19[:,:],norm=LogNorm(),levels=levels_mass*10, cmap='Oranges', alpha=0.75)
ax10.set_yscale('log')
ax10.set_ylim(100.,min(p)/1.e6)
ax10.set_xlim(min(longitude),max(longitude))
ax10.text(100,1e-5,'1900 K')

im2 = ax11.contourf(X2,Y2,al_mass_20[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Greens', alpha=0.75)
im1 = ax11.contourf(X1,Y1,ti_mass_20[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Blues', alpha=0.75)
#im4 = ax11.contourf(X4,Y4,mg_mass_20[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
im3 = ax11.contourf(X3,Y3,fe_mass_20[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Reds', alpha=0.75)
im5 = ax11.contourf(X5,Y5,cr_mass_20[:,:],norm=LogNorm(),levels=levels_mass*10, cmap='Oranges', alpha=0.75)
ax11.set_yscale('log')
ax11.set_ylim(100.,min(p)/1.e6)
ax11.set_xlim(min(longitude),max(longitude))
ax11.text(100,1e-5,'2000 K')

im2 = ax12.contourf(X2,Y2,al_mass_21[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Greens', alpha=0.75)
im1 = ax12.contourf(X1,Y1,ti_mass_21[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Blues', alpha=0.75)
#im4 = ax12.contourf(X4,Y4,mg_mass_21[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
im3 = ax12.contourf(X3,Y3,fe_mass_21[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Reds', alpha=0.75)
im5 = ax12.contourf(X5,Y5,cr_mass_21[:,:],norm=LogNorm(),levels=levels_mass*10, cmap='Oranges', alpha=0.75)
ax12.set_yscale('log')
ax12.set_ylim(100.,min(p)/1.e6)
ax12.set_xlim(min(longitude),max(longitude))
ax12.text(100,1e-5,'2100 K')

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.83, 0.15, 0.02, 0.4])
cbar = fig.colorbar(im1, cax=cbar_ax)
cbar.set_label(r'Mass Density [g cm$^{-3}$]', fontsize=fontsized, rotation=-90, labelpad=30)

plt.subplots_adjust(wspace=0.02, hspace=0.035)
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.xlabel(r"Longitude [degree]")
plt.ylabel(r"Pressure [bar]", labelpad = 25)

ax1.tick_params(axis='both', which='major', labelsize=fontsized)
ax2.tick_params(axis='both', which='major', labelsize=fontsized)
ax3.tick_params(axis='both', which='major', labelsize=fontsized)
ax4.tick_params(axis='both', which='major', labelsize=fontsized)
ax5.tick_params(axis='both', which='major', labelsize=fontsized)
ax6.tick_params(axis='both', which='major', labelsize=fontsized)
ax7.tick_params(axis='both', which='major', labelsize=fontsized)
ax8.tick_params(axis='both', which='major', labelsize=fontsized)
ax9.tick_params(axis='both', which='major', labelsize=fontsized)
ax10.tick_params(axis='both', which='major', labelsize=fontsized)
ax11.tick_params(axis='both', which='major', labelsize=fontsized)
ax12.tick_params(axis='both', which='major', labelsize=fontsized)

#plt.savefig(path+'/layers_mass_all.jpg', format = 'jpg', bbox_inches='tight')
plt.savefig(path+'/layers_mass_all.pdf', format = 'pdf', bbox_inches='tight')
plt.close()

fig = plt.figure()
fig, ((ax1,ax2),(ax3,ax4),(ax5,ax6),(ax7,ax8),(ax9,ax10),(ax11,ax12)) = plt.subplots(nrows = 6, ncols = 2,sharey=True, sharex= True, figsize=(15,32))

levels_mass = np.array([1e-12,1e-10,1e-8,1e-6,1e-4,1e-2])/1e7


im1 = ax1.contourf(X4,Y4,mg_mass_10[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
#ax1.plot(longitude,mg_cond_temp_10, c='k', lw=2)
ax1.set_yscale('log')
ax1.set_ylim(100.,min(p)/1.e6)
ax1.set_xlim(min(longitude),max(longitude))
ax1.text(100,1e-5,'1000 K')


im1 = ax2.contourf(X4,Y4,mg_mass_11[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
#ax2.plot(longitude,mg_cond_temp_11, c='k', lw=2)
ax2.set_yscale('log')
ax2.set_ylim(100.,min(p)/1.e6)
ax2.set_xlim(min(longitude),max(longitude))
ax2.text(100,1e-5,'1100 K')

im1 = ax3.contourf(X2,Y2,mg_mass_12[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
#ax3.plot(longitude,mg_cond_temp_12, c='k', lw=2)
ax3.set_yscale('log')
ax3.set_ylim(100.,min(p)/1.e6)
ax3.set_xlim(min(longitude),max(longitude))
ax3.text(100,1e-5,'1200 K')

im1 = ax4.contourf(X2,Y2,mg_mass_13[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
#ax4.plot(longitude,mg_cond_temp_13, c='k', lw=2)
ax4.set_yscale('log')
ax4.set_ylim(100.,min(p)/1.e6)
ax4.set_xlim(min(longitude),max(longitude))
ax4.text(100,1e-5,'1300 K')

im1 = ax5.contourf(X4,Y4,mg_mass_14[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
#ax5.plot(longitude,mg_cond_temp_14, c='k', lw=2)
ax5.set_yscale('log')
ax5.set_ylim(100.,min(p)/1.e6)
ax5.set_xlim(min(longitude),max(longitude))
ax5.text(100,1e-5,'1400 K')

im1 = ax6.contourf(X4,Y4,mg_mass_15[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
#ax6.plot(longitude,mg_cond_temp_15, c='k', lw=2)
ax6.set_yscale('log')
ax6.set_ylim(100.,min(p)/1.e6)
ax6.set_xlim(min(longitude),max(longitude))
ax6.text(100,1e-5,'1500 K')

im1 = ax7.contourf(X4,Y4,mg_mass_16[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
#ax7.plot(longitude,mg_cond_temp_16, c='k', lw=2)
ax7.set_yscale('log')
ax7.set_ylim(100.,min(p)/1.e6)
ax7.set_xlim(min(longitude),max(longitude))
ax7.text(100,1e-5,'1600 K')

im1 = ax8.contourf(X4,Y4,mg_mass_17[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
#ax8.plot(longitude,mg_cond_temp_17, c='k', lw=2)
ax8.set_yscale('log')
ax8.set_ylim(100.,min(p)/1.e6)
ax8.set_xlim(min(longitude),max(longitude))
ax8.text(100,1e-5,'1700 K')

im1 = ax9.contourf(X4,Y4,mg_mass_18[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
#ax9.plot(longitude,mg_cond_temp_18, c='k', lw=2)
ax9.set_yscale('log')
ax9.set_ylim(100.,min(p)/1.e6)
ax9.set_xlim(min(longitude),max(longitude))
ax9.text(100,1e-5,'1800 K')

im1 = ax10.contourf(X4,Y4,mg_mass_19[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
#ax10.plot(longitude,mg_cond_temp_19, c='k', lw=2)
ax10.set_yscale('log')
ax10.set_ylim(100.,min(p)/1.e6)
ax10.set_xlim(min(longitude),max(longitude))
ax10.text(100,1e-5,'1900 K')

im1 = ax11.contourf(X4,Y4,mg_mass_20[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
#ax11.plot(longitude,mg_cond_temp_20, c='k', lw=2)
ax11.set_yscale('log')
ax11.set_ylim(100.,min(p)/1.e6)
ax11.set_xlim(min(longitude),max(longitude))
ax11.text(100,1e-5,'2000 K')

im1 = ax12.contourf(X4,Y4,mg_mass_21[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
#ax12.plot(longitude,mg_cond_temp_21, c='k', lw=2)
ax12.set_yscale('log')
ax12.set_ylim(100.,min(p)/1.e6)
ax12.set_xlim(min(longitude),max(longitude))
ax12.text(100,1e-5,'2100 K')


fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.83, 0.15, 0.02, 0.4])
cbar = fig.colorbar(im1, cax=cbar_ax)
cbar.set_label(r'Mass Density [g cm$^{-3}$]', fontsize=fontsized, rotation=-90, labelpad=30)


plt.subplots_adjust(wspace=0.02, hspace=0.035)
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.xlabel(r"Longitude [degree]")
plt.ylabel(r"Pressure [bar]", labelpad = 25)

#plt.savefig(path+'/layers_mass_all_mg.jpg', format = 'jpg', bbox_inches='tight')
plt.savefig(path+'/layers_mass_all_mg.pdf', format = 'pdf', bbox_inches='tight')
plt.close()

