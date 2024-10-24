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
name = args[1]
path = args[2]
tp_file_east = '/data/groups/zhang/dkpowell/working_dir/2D_setup_zen/Lux_2D_setup/PT_grid_winds/exotransmit_inputs/exotransmit_input_'+str(name)+'east.txt'
tp_file_west = '/data/groups/zhang/dkpowell/working_dir/2D_setup_zen/Lux_2D_setup/PT_grid_winds/exotransmit_inputs/exotransmit_input_'+str(name)+'west.txt'

# Check if folder where the plots are generated exists; if 
# not, make one. 
if not os.path.isdir(path):
	os.makedirs(path)
	
# Read in given output file and plot results. 
print('Read in output values file')
infile=open(file,'r')
line = infile.readline().split()

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
			

al_east = np.loadtxt('east_cloud_properties'+str(name)+'_al.txt')
al_east = al_east.T #nz,nbin
cr_east = np.loadtxt('east_cloud_properties'+str(name)+'_cr.txt')
cr_east = cr_east.T #nz,nbin
fe_east = np.loadtxt('east_cloud_properties'+str(name)+'_fe.txt')
fe_east = fe_east.T #nz,nbin
kc_east = np.loadtxt('east_cloud_properties'+str(name)+'_kc.txt')
kc_east = kc_east.T #nz,nbin
mg_east = np.loadtxt('east_cloud_properties'+str(name)+'_mg.txt')
mg_east = mg_east.T #nz,nbin
mn_east = np.loadtxt('east_cloud_properties'+str(name)+'_mn.txt')
mn_east = mn_east.T #nz,nbin
na_east = np.loadtxt('east_cloud_properties'+str(name)+'_na.txt')
na_east = na_east.T #nz,nbin
ti_east = np.loadtxt('east_cloud_properties'+str(name)+'_ti.txt')
ti_east = ti_east.T #nz,nbin
zn_east = np.loadtxt('east_cloud_properties'+str(name)+'_zn.txt')
zn_east = zn_east.T #nz,nbin

al_west = np.loadtxt('west_cloud_properties'+str(name)+'_al.txt')
al_west = al_west.T #nz,nbin
cr_west = np.loadtxt('west_cloud_properties'+str(name)+'_cr.txt')
cr_west = cr_west.T #nz,nbin
fe_west = np.loadtxt('west_cloud_properties'+str(name)+'_fe.txt')
fe_west = fe_west.T #nz,nbin
kc_west = np.loadtxt('west_cloud_properties'+str(name)+'_kc.txt')
kc_west = kc_west.T #nz,nbin
mg_west = np.loadtxt('west_cloud_properties'+str(name)+'_mg.txt')
mg_west = mg_west.T #nz,nbin
mn_west = np.loadtxt('west_cloud_properties'+str(name)+'_mn.txt')
mn_west = mn_west.T #nz,nbin
na_west = np.loadtxt('west_cloud_properties'+str(name)+'_na.txt')
na_west = na_west.T #nz,nbin
ti_west = np.loadtxt('west_cloud_properties'+str(name)+'_ti.txt')
ti_west = ti_west.T #nz,nbin
zn_west = np.loadtxt('west_cloud_properties'+str(name)+'_zn.txt')
zn_west = zn_west.T #nz,nbin


wavelength,extinct_ti = np.loadtxt('rf_data/ti_rad_1.txt', unpack=True)


od_ti_east  = zeros((nz,nbin,len(wavelength)))
od_mg_east  = zeros((nz,nbin,len(wavelength)))
od_fe_east  = zeros((nz,nbin,len(wavelength)))
od_cr_east  = zeros((nz,nbin,len(wavelength)))
od_mn_east  = zeros((nz,nbin,len(wavelength)))
od_na_east  = zeros((nz,nbin,len(wavelength)))
od_zn_east  = zeros((nz,nbin,len(wavelength)))
od_k_east = zeros((nz,nbin,len(wavelength)))
od_al_east  = zeros((nz,nbin,len(wavelength)))

od_tot_east = zeros((nz,nbin,len(wavelength)))

od_ti_west  = zeros((nz,nbin,len(wavelength)))
od_mg_west  = zeros((nz,nbin,len(wavelength)))
od_fe_west  = zeros((nz,nbin,len(wavelength)))
od_cr_west  = zeros((nz,nbin,len(wavelength)))
od_mn_west  = zeros((nz,nbin,len(wavelength)))
od_na_west  = zeros((nz,nbin,len(wavelength)))
od_zn_west  = zeros((nz,nbin,len(wavelength)))
od_k_west = zeros((nz,nbin,len(wavelength)))
od_al_west  = zeros((nz,nbin,len(wavelength)))

od_tot_west = zeros((nz,nbin,len(wavelength)))


sum_od_ti = zeros((nz,len(wavelength)))
sum_od_mg = zeros((nz,len(wavelength)))
sum_od_fe = zeros((nz,len(wavelength)))
sum_od_cr = zeros((nz,len(wavelength)))
sum_od_mn = zeros((nz,len(wavelength)))
sum_od_na = zeros((nz,len(wavelength)))
sum_od_zn = zeros((nz,len(wavelength)))
sum_od_k = zeros((nz,len(wavelength)))
sum_od_al = zeros((nz,len(wavelength)))

tot_cloud_opaque = zeros((nz,len(wavelength)))

sum_od_ti_east = zeros((nz,len(wavelength)))
sum_od_mg_east = zeros((nz,len(wavelength)))
sum_od_fe_east = zeros((nz,len(wavelength)))
sum_od_cr_east = zeros((nz,len(wavelength)))
sum_od_mn_east = zeros((nz,len(wavelength)))
sum_od_na_east = zeros((nz,len(wavelength)))
sum_od_zn_east = zeros((nz,len(wavelength)))
sum_od_k_east= zeros((nz,len(wavelength)))
sum_od_al_east = zeros((nz,len(wavelength)))

tot_cloud_opaque_east = zeros((nz,len(wavelength)))


sum_od_ti_west = zeros((nz,len(wavelength)))
sum_od_mg_west = zeros((nz,len(wavelength)))
sum_od_fe_west = zeros((nz,len(wavelength)))
sum_od_cr_west = zeros((nz,len(wavelength)))
sum_od_mn_west = zeros((nz,len(wavelength)))
sum_od_na_west = zeros((nz,len(wavelength)))
sum_od_zn_west = zeros((nz,len(wavelength)))
sum_od_k_west= zeros((nz,len(wavelength)))
sum_od_al_west = zeros((nz,len(wavelength)))

tot_cloud_opaque_west = zeros((nz,len(wavelength)))

	
ssa_tot = zeros((nz,nbin,len(wavelength)))

g_tot = zeros((nz,nbin,len(wavelength)))

extinct_ti = zeros((len(wavelength),nbin))
extinct_mg = zeros((len(wavelength),nbin))
extinct_fe = zeros((len(wavelength),nbin))
extinct_cr = zeros((len(wavelength),nbin))
extinct_mn = zeros((len(wavelength),nbin))
extinct_na = zeros((len(wavelength),nbin))
extinct_zn = zeros((len(wavelength),nbin))
extinct_k = zeros((len(wavelength),nbin))
extinct_al = zeros((len(wavelength),nbin))

ssa_ti = zeros((len(wavelength),nbin))
ssa_mg = zeros((len(wavelength),nbin))
ssa_fe = zeros((len(wavelength),nbin))
ssa_cr = zeros((len(wavelength),nbin))
ssa_mn = zeros((len(wavelength),nbin))
ssa_na = zeros((len(wavelength),nbin))
ssa_zn = zeros((len(wavelength),nbin))
ssa_k = zeros((len(wavelength),nbin))
ssa_al = zeros((len(wavelength),nbin))

g_ti = zeros((len(wavelength),nbin))
g_mg = zeros((len(wavelength),nbin))
g_fe = zeros((len(wavelength),nbin))
g_cr = zeros((len(wavelength),nbin))
g_mn = zeros((len(wavelength),nbin))
g_na = zeros((len(wavelength),nbin))
g_zn = zeros((len(wavelength),nbin))
g_k = zeros((len(wavelength),nbin))
g_al = zeros((len(wavelength),nbin))

iwave,ibin,ext_tio2, tio2_ssa, tio2_g =    np.loadtxt('miex_data_morebins/tio2_opticalproperties.txt', unpack = True)
iwave,ibin,ext_mns, mns_ssa, mns_g =       np.loadtxt('miex_data_morebins/mns_opticalproperties.txt', unpack = True)
iwave,ibin,ext_cr, cr_ssa, cr_g =          np.loadtxt('miex_data_morebins/cr_opticalproperties.txt', unpack = True)
iwave,ibin,ext_al2o3, al2o3_ssa, al2o3_g = np.loadtxt('miex_data_morebins/al2o3_opticalproperties.txt', unpack = True)
iwave,ibin,ext_na2s, na2s_ssa, na2s_g =    np.loadtxt('miex_data_morebins/na2s_opticalproperties.txt', unpack = True)
iwave,ibin,ext_mg2sio4, mg2sio4_ssa, mg2sio4_g = np.loadtxt('miex_data_morebins/mg2sio4_opticalproperties.txt', unpack = True)
iwave,ibin,ext_kcl, kcl_ssa, kcl_g =       np.loadtxt('miex_data_morebins/kcl_opticalproperties.txt', unpack = True)
iwave,ibin,ext_zns, zns_ssa, zns_g =       np.loadtxt('miex_data_morebins/zns_opticalproperties.txt', unpack = True)
iwave,ibin,ext_fe, fe_ssa, fe_g =          np.loadtxt('miex_data_morebins/fe_opticalproperties.txt', unpack = True)

count = 0

for idw,w in enumerate(wavelength):
	for j in range(nbin):
		extinct_ti[idw,j] = np.abs(ext_tio2[count])
		extinct_mg[idw,j] = np.abs(ext_mg2sio4[count])
		extinct_mn[idw,j] = np.abs(ext_mns[count])
		extinct_na[idw,j] = np.abs(ext_na2s[count])
		extinct_k[idw,j] = np.abs(ext_kcl[count])
		extinct_fe[idw,j] = np.abs(ext_fe[count])
		extinct_cr[idw,j] = np.abs(ext_cr[count])
		extinct_al[idw,j] = np.abs(ext_al2o3[count])
		extinct_zn[idw,j] = np.abs(ext_zns[count])
		
		ssa_ti[idw,j] = np.abs(tio2_ssa[count])
		ssa_mg[idw,j] = np.abs(mg2sio4_ssa[count])
		ssa_mn[idw,j] = np.abs(mns_ssa[count])
		ssa_na[idw,j] = np.abs(na2s_ssa[count])
		ssa_k[idw,j]  = np.abs(kcl_ssa[count])
		ssa_fe[idw,j] = np.abs(fe_ssa[count])
		ssa_cr[idw,j] = np.abs(cr_ssa[count])
		ssa_al[idw,j] = np.abs(al2o3_ssa[count])
		ssa_zn[idw,j] = np.abs(zns_ssa[count])
		
		g_ti[idw,j] = tio2_g[count]
		g_mg[idw,j] = mg2sio4_g[count]
		g_mn[idw,j] = mns_g[count]
		g_na[idw,j] = na2s_g[count]
		g_k[idw,j]  = kcl_g[count]
		g_fe[idw,j] = fe_g[count]
		g_cr[idw,j] = cr_g[count]
		g_al[idw,j] = al2o3_g[count]
		g_zn[idw,j] = zns_g[count]
		count +=1
		

##########################################################################################
##########################################################################################	
print('BIG LOOP TIME!')
##########################################################################################
########################################################################################## 

for idw,w in enumerate(wavelength):
	for k in range(nz): 
		for j in range(nbin):	
			od_ti_east[k,j,idw] = ti_east[k,j]*extinct_ti[idw,j]
			od_mg_east[k,j,idw] = mg_east[k,j]*extinct_mg[idw,j]
			od_mn_east[k,j,idw] = mn_east[k,j]*extinct_mn[idw,j]
			od_na_east[k,j,idw] = na_east[k,j]*extinct_na[idw,j]
			od_zn_east[k,j,idw] = zn_east[k,j]*extinct_zn[idw,j]
			od_k_east[k,j,idw] = kc_east[k,j]*extinct_k[idw,j]
			od_al_east[k,j,idw] = al_east[k,j]*extinct_al[idw,j]
			od_fe_east[k,j,idw] = fe_east[k,j]*extinct_fe[idw,j]
			od_cr_east[k,j,idw] = cr_east[k,j]*extinct_cr[idw,j]
			
			od_tot_east[k,j,idw] = od_ti_east[k,j,idw]+od_mg_east[k,j,idw]+od_mn_east[k,j,idw]+od_na_east[k,j,idw]+od_zn_east[k,j,idw]+od_k_east[k,j,idw]+od_al_east[k,j,idw]+od_fe_east[k,j,idw]+od_cr_east[k,j,idw]
			
			od_ti_west[k,j,idw] = ti_west[k,j]*extinct_ti[idw,j]
			od_mg_west[k,j,idw] = mg_west[k,j]*extinct_mg[idw,j]
			od_mn_west[k,j,idw] = mn_west[k,j]*extinct_mn[idw,j]
			od_na_west[k,j,idw] = na_west[k,j]*extinct_na[idw,j]
			od_zn_west[k,j,idw] = zn_west[k,j]*extinct_zn[idw,j]
			od_k_west[k,j,idw] = kc_west[k,j]*extinct_k[idw,j]
			od_al_west[k,j,idw] = al_west[k,j]*extinct_al[idw,j]
			od_fe_west[k,j,idw] = fe_west[k,j]*extinct_fe[idw,j]
			od_cr_west[k,j,idw] = cr_west[k,j]*extinct_cr[idw,j]
			
			od_tot_west[k,j,idw] = od_ti_west[k,j,idw]+od_mg_west[k,j,idw]+od_mn_west[k,j,idw]+od_na_west[k,j,idw]+od_zn_west[k,j,idw]+od_k_west[k,j,idw]+od_al_west[k,j,idw]+od_fe_west[k,j,idw]+od_cr_west[k,j,idw]
			
			
		
		sum_od_ti_east[k,idw] = sum(od_ti_east[k,:,idw])	
		sum_od_mg_east[k,idw] = sum(od_mg_east[k,:,idw])	
		sum_od_mn_east[k,idw] = sum(od_mn_east[k,:,idw])	
		sum_od_na_east[k,idw] = sum(od_na_east[k,:,idw])	
		sum_od_zn_east[k,idw] = sum(od_zn_east[k,:,idw])	
		sum_od_k_east[k,idw] = sum(od_k_east[k,:,idw])	
		sum_od_al_east[k,idw] = sum(od_al_east[k,:,idw])	
		sum_od_fe_east[k,idw] = sum(od_fe_east[k,:,idw])	
		sum_od_cr_east[k,idw] = sum(od_cr_east[k,:,idw])	
		
		
		
		tot_cloud_opaque_east[k,idw] = sum_od_ti_east[k,idw]+sum_od_mg_east[k,idw]+sum_od_mn_east[k,idw]+sum_od_na_east[k,idw]+sum_od_zn_east[k,idw]+sum_od_k_east[k,idw]+sum_od_al_east[k,idw]+sum_od_fe_east[k,idw]+sum_od_cr_east[k,idw]

		sum_od_ti_west[k,idw] = sum(od_ti_west[k,:,idw])	
		sum_od_mg_west[k,idw] = sum(od_mg_west[k,:,idw])	
		sum_od_mn_west[k,idw] = sum(od_mn_west[k,:,idw])	
		sum_od_na_west[k,idw] = sum(od_na_west[k,:,idw])	
		sum_od_zn_west[k,idw] = sum(od_zn_west[k,:,idw])	
		sum_od_k_west[k,idw] = sum(od_k_west[k,:,idw])	
		sum_od_al_west[k,idw] = sum(od_al_west[k,:,idw])	
		sum_od_fe_west[k,idw] = sum(od_fe_west[k,:,idw])	
		sum_od_cr_west[k,idw] = sum(od_cr_west[k,:,idw])	
		
		
		
		tot_cloud_opaque_west[k,idw] = sum_od_ti_west[k,idw]+sum_od_mg_west[k,idw]+sum_od_mn_west[k,idw]+sum_od_na_west[k,idw]+sum_od_zn_west[k,idw]+sum_od_k_west[k,idw]+sum_od_al_west[k,idw]+sum_od_fe_west[k,idw]+sum_od_cr_west[k,idw]


f_tot_opaque_east = interp2d(wavelength,p/1.e6,tot_cloud_opaque_east)

f_tot_opaque_west = interp2d(wavelength,p/1.e6,tot_cloud_opaque_west)


i,pressure_east,temperature_east = np.loadtxt(tp_file_east, unpack=True, skiprows = 1)
pressure_east*=1e-5

i,pressure_west,temperature_west = np.loadtxt(tp_file_west, unpack=True, skiprows = 1)
pressure_west*=1e-5

print(p/1.e6, pressure_east)

output_east = f_tot_opaque_east(wavelength,pressure_east)
output_east *= 100. #opacity is in units of m^-1

output_west = f_tot_opaque_west(wavelength,pressure_west)
output_west *= 100. #opacity is in units of m^-1


np.savetxt(path+'/cloud_opacities_east_'+name+'.txt', output_east.T, fmt='%.3e')

np.savetxt(path+'/cloud_opacities_west_'+name+'.txt', output_west.T, fmt='%.3e')


		
			