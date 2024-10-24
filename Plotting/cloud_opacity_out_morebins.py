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

#./cloud_opacity_out_morebins.py hjt_west_1800.txt ../T_P/t_p_1800_west_limb.dat 1800_west opacities/ > log.txt &

#./cloud_opacity_out_morebins.py ../../../New_ExoCARMA/apr7_half_latentheat/run/carma/hjt_west_2000.txt ../T_P/t_p_2000_west_limb.dat 2000_west apr11_halflatent_opacities/ > log.txt &
#./cloud_opacity_out_morebins.py ../../../New_ExoCARMA/apr7_half_latentheat/run/carma/hjt_east_2000.txt ../T_P/t_p_2000_east_limb.dat 2000_east apr11_halflatent_opacities/ > log.txt &
#./cloud_opacity_out_morebins.py ../../../New_ExoCARMA/apr7_half_latentheat/run/carma/hjt_pole_2000.txt ../T_P/t_p_2000_pole.dat 2000_pole apr11_halflatent_opacities/ > log.txt &

# Read in inputs from command line
opts, args = getopt.getopt(sys.argv[1:], "h")
file = args[0]
tp_file = args[1]
name = args[2]
path = args[3]

# Check if folder where the plots are generated exists; if 
# not, make one. 
if not os.path.isdir(path):
	os.makedirs(path)
	
# Read in given output file and plot results. 
print('Read in output file')
data=open(file,'r')
line = data.readline().split()

nz,ngroup,nelem,nbin,ngas,nstep,iskip = map(int,line)
r,ms,dr,rl,ru = zeros((nbin,ngroup)),zeros((nbin,ngroup)),zeros((nbin,ngroup)),zeros((nbin,ngroup)),zeros((nbin,ngroup))
p,t,z = zeros(nz),zeros(nz),zeros(nz)


ntime = int(nstep/iskip)

if name == 'west_big':
	ntime = int(789000/iskip)
if name == 'day_big':
	ntime = int(932000/iskip)

time = zeros(ntime)

tio2, het_fe, het_mg2sio4, pure_fe, het_cr, cr, mns, na2s,zns,kcl,al2o3 = zeros((nz,nbin,ntime)),zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime))
mass_tio2, mass_het_fe, mass_het_mg2sio4, mass_pure_fe, mass_het_cr, mass_cr, mass_mns, mass_na2s, mass_zns, mass_kcl, mass_al2o3 = zeros((nz,nbin,ntime)),zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime))

t_tio2, t_het_fe, t_het_mg2sio4, t_pure_fe, t_het_cr, t_cr, t_mns, t_na2s,t_zns,t_kcl,t_al2o3 = zeros((nz,nbin)),zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin))
t_mass_tio2, t_mass_het_fe, t_mass_het_mg2sio4, t_mass_pure_fe, t_mass_het_cr, t_mass_cr, t_mass_mns, t_mass_na2s, t_mass_zns, t_mass_kcl, t_mass_al2o3 = zeros((nz,nbin)),zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin))

p_ti, p_fe, p_mg, p_cr, p_mn, p_na,p_zn,p_k,p_al = zeros((nz,nbin)),zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin))
pm_ti, pm_fe, pm_mg, pm_cr, pm_mn, pm_na,pm_zn,pm_k,pm_al = zeros((nz,nbin)),zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin))
p_pure_cr,p_pure_fe = zeros((nz,nbin)),zeros((nz,nbin))
pm_pure_cr,pm_pure_fe = zeros((nz,nbin)),zeros((nz,nbin))

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
		line = data.readline().split()
		r[j,i] = float(line[2])
		ms[j,i] = float(line[3])
		dr[j,i] = float(line[4])
		rl[j,i] = float(line[5])
		ru[j,i] = float(line[6])
for i in range(nz):
	line = data.readline().split()
	z[i] = float(line[1])
	p[i] = float(line[3])
	t[i] = float(line[4])
			
np.savetxt('new_ti_rad.txt', r[:,0].T)
np.savetxt('new_al_rad.txt', r[:,1].T)
np.savetxt('new_fe_rad.txt', r[:,2].T)
np.savetxt('new_mg_rad.txt', r[:,3].T)
np.savetxt('new_cr_rad.txt', r[:,4].T)
np.savetxt('new_mn_rad.txt', r[:,5].T)
np.savetxt('new_na_rad.txt', r[:,6].T)
np.savetxt('new_pfe_rad.txt', r[:,7].T)
np.savetxt('new_pcr_rad.txt', r[:,8].T)
np.savetxt('new_k_rad.txt', r[:,9].T)
np.savetxt('new_zn_rad.txt', r[:,10].T)

for i in range(ntime-1):
	time[i] = float(data.readline())
	for j in range(nbin):
		for k in range(nz):
			line = data.readline().split()
			
			tio2[k,j,i] = float(line[2])/(math.log(ru[j,0]) - math.log(rl[j,0]))
			mass_tio2[k,j,i] = tio2[k,j,i]*(M_TIO2*((4./3.)*np.pi*(r[j,0]*1e-4)**3.))
			
			al2o3[k,j,i] = float(line[3])/(math.log(ru[j,1]) - math.log(rl[j,1]))
			mass_al2o3[k,j,i] = al2o3[k,j,i]*(M_AL2O3*((4./3.)*np.pi*(r[j,1]*1e-4)**3.))
			
			het_fe[k,j,i] = float(line[4])/(math.log(ru[j,2]) - math.log(rl[j,2]))
			mass_het_fe[k,j,i] = het_fe[k,j,i]*(M_FE*((4./3.)*np.pi*(r[j,2]*1e-4)**3.))
			
			het_mg2sio4[k,j,i] = float(line[5])/(math.log(ru[j,3]) - math.log(rl[j,3]))
			mass_het_mg2sio4[k,j,i] = het_mg2sio4[k,j,i]*(M_MG2SIO4*((4./3.)*np.pi*(r[j,3]*1e-4)**3.))
			
			het_cr[k,j,i] = float(line[6])/(math.log(ru[j,4]) - math.log(rl[j,4]))
			mass_het_cr[k,j,i] = het_cr[k,j,i]*(M_CR*((4./3.)*np.pi*(r[j,4]*1e-4)**3.))
			
			mns[k,j,i] = float(line[7])/(math.log(ru[j,5]) - math.log(rl[j,5]))
			mass_mns[k,j,i] = mns[k,j,i]*(M_MNS*((4./3.)*np.pi*(r[j,5]*1e-4)**3.))
			
			na2s[k,j,i] = float(line[8])/(math.log(ru[j,6]) - math.log(rl[j,6]))
			mass_na2s[k,j,i] = na2s[k,j,i]*(M_NA2S*((4./3.)*np.pi*(r[j,6]*1e-4)**3.))
			
			pure_fe[k,j,i] = float(line[9])/(math.log(ru[j,7]) - math.log(rl[j,7]))
			mass_pure_fe[k,j,i] = pure_fe[k,j,i]*(M_FE*((4./3.)*np.pi*(r[j,7]*1e-4)**3.))
			
			cr[k,j,i] = float(line[10])/(math.log(ru[j,8]) - math.log(rl[j,8]))
			mass_cr[k,j,i] = cr[k,j,i]*(M_CR*((4./3.)*np.pi*(r[j,8]*1e-4)**3.))
			
			kcl[k,j,i] = float(line[11])/(math.log(ru[j,9]) - math.log(rl[j,9]))
			mass_kcl[k,j,i] = kcl[k,j,i]*(M_KCL*((4./3.)*np.pi*(r[j,9]*1e-4)**3.))
			
			zns[k,j,i] = float(line[12])/(math.log(ru[j,10]) - math.log(rl[j,10]))
			mass_zns[k,j,i] = zns[k,j,i]*(M_ZNS*((4./3.)*np.pi*(r[j,10]*1e-4)**3.))
					
	if i > ntime-10 :
		count += 1
			
		t_tio2[:,:] += tio2[:,:,i]
		t_mass_tio2[:,:] += mass_tio2[:,:,i]
			
		t_het_fe[:,:] += het_fe[:,:,i]
		t_mass_het_fe[:,:] += mass_het_fe[:,:,i]
			
		t_het_mg2sio4[:,:] += het_mg2sio4[:,:,i]
		t_mass_het_mg2sio4[:,:] += mass_het_mg2sio4[:,:,i]
			
		t_pure_fe[:,:] += pure_fe[:,:,i]
		t_mass_pure_fe[:,:] += mass_pure_fe[:,:,i]
			
		t_het_cr[:,:] += het_cr[:,:,i]
		t_mass_het_cr[:,:] += mass_het_cr[:,:,i]
			
		t_cr[:,:] += cr[:,:,i]
		t_mass_cr[:,:] += mass_cr[:,:,i]
			
		t_mns[:,:] += mns[:,:,i]
		t_mass_mns[:,:] += mass_mns[:,:,i]
			
		t_na2s[:,:] += na2s[:,:,i]
		t_mass_na2s[:,:] += mass_na2s[:,:,i]
			
		t_zns[:,:] += zns[:,:,i]
		t_mass_zns[:,:] += mass_zns[:,:,i]
			
		t_kcl[:,:] += kcl[:,:,i]
		t_mass_kcl[:,:] += mass_kcl[:,:,i]
			
		t_al2o3[:,:] += al2o3[:,:,i]
		t_mass_al2o3[:,:] += mass_al2o3[:,:,i]
						
t_tio2 /= count 
t_mass_tio2 /= count 
t_het_fe /= count
t_mass_het_fe /= count
t_het_mg2sio4 /= count
t_mass_het_mg2sio4 /= count
t_pure_fe /= count
t_mass_pure_fe /= count
t_het_cr /= count
t_mass_het_cr /= count
t_cr /= count
t_mass_cr /= count
t_mns /= count
t_mass_mns /= count
t_na2s /= count
t_mass_na2s /= count
t_zns /= count
t_mass_zns /= count
t_kcl /= count
t_mass_kcl /= count
t_al2o3 /= count
t_mass_al2o3 /= count

#print np.max(al2o3), np.max(het_mg2sio4), np.max(het_fe+pure_fe), np.max(het_cr+cr), np.max(tio2), np.max(mns)
#print np.max(mass_al2o3), np.max(mass_het_mg2sio4), np.max(mass_het_fe+mass_pure_fe), np.max(mass_het_cr+mass_cr), np.max(mass_tio2), np.max(mass_mns)

num_cut = 1e-5
mass_cut = 1e-15

for j in range(nbin):
	for k in range(nz):
		if t_tio2[k,j] > num_cut:
			p_ti[k,j] = t_tio2[k,j]
		else:
			p_ti[k,j] = 0.
			
		if t_het_fe[k,j]+t_pure_fe[k,j] > num_cut:
			p_fe[k,j] = t_het_fe[k,j]
			p_pure_fe[k,j] = t_pure_fe[k,j]
		else:
			p_fe[k,j] = 0.
			p_pure_fe[k,j] = 0.
		
		if t_het_cr[k,j]+t_cr[k,j] > num_cut:
			p_cr[k,j] = t_het_cr[k,j]
			p_pure_cr[k,j] = t_cr[k,j]
		else:
			p_cr[k,j] = 0.
			p_pure_cr[k,j] = 0.
		
		if t_het_mg2sio4[k,j] > num_cut:
			p_mg[k,j] = t_het_mg2sio4[k,j]
		else:
			p_mg[k,j] = 0.
	
		if t_mns[k,j] > num_cut:
			p_mn[k,j] = t_mns[k,j]
		else:
			p_mn[k,j] = 0.
			
		if t_na2s[k,j] > num_cut:
			p_na[k,j] = t_na2s[k,j]
		else:
			p_na[k,j] = 0.
			
		if t_zns[k,j] > num_cut:
			p_zn[k,j] = t_zns[k,j]
		else:
			p_zn[k,j] = 0.
		
		if t_kcl[k,j] > num_cut:
			p_k[k,j] = t_kcl[k,j]
		else:
			p_k[k,j] = 0.
			
		if t_al2o3[k,j] > num_cut:
			p_al[k,j] = t_al2o3[k,j]
		else:
			p_al[k,j] = 0.
		
		if t_mass_tio2[k,j] > mass_cut:
			pm_ti[k,j] = t_mass_tio2[k,j]
		else:
			pm_ti[k,j] = 0.
			
		if t_mass_het_fe[k,j]+t_mass_pure_fe[k,j] > mass_cut:
			pm_fe[k,j] = t_mass_het_fe[k,j]
			pm_pure_fe[k,j] = t_mass_pure_fe[k,j]
		else:
			pm_fe[k,j] = 0.
			pm_pure_fe[k,j] = 0.
		
		if t_mass_het_cr[k,j]+t_mass_cr[k,j] > mass_cut:
			pm_cr[k,j] = t_mass_het_cr[k,j]
			pm_pure_cr[k,j] = t_mass_cr[k,j]
		else:
			pm_cr[k,j] = 0.
			pm_pure_cr[k,j] = 0.
		
		if t_mass_het_mg2sio4[k,j] > mass_cut:
			pm_mg[k,j] = t_mass_het_mg2sio4[k,j]
		else:
			pm_mg[k,j] = 0.
	
		if t_mass_mns[k,j] > mass_cut:
			pm_mn[k,j] = t_mass_mns[k,j]
		else:
			pm_mn[k,j] = 0.
			
		if t_mass_na2s[k,j] > mass_cut:
			pm_na[k,j] = t_mass_na2s[k,j]
		else:
			pm_na[k,j] = 0.
			
		if t_mass_zns[k,j] > mass_cut:
			pm_zn[k,j] = t_mass_zns[k,j]
		else:
			pm_zn[k,j] = 0.
		
		if t_mass_kcl[k,j] > mass_cut:
			pm_k[k,j] = t_mass_kcl[k,j]
		else:
			pm_k[k,j] = 0.
			
		if t_mass_al2o3[k,j] > mass_cut:
			pm_al[k,j] = t_mass_al2o3[k,j]
		else:
			pm_al[k,j] = 0.
			
wavelength,extinct_ti = np.loadtxt('rf_data/ti_rad_1.txt', unpack=True)


od_ti = zeros((nz,nbin,len(wavelength)))
od_mg = zeros((nz,nbin,len(wavelength)))
od_fe = zeros((nz,nbin,len(wavelength)))
od_cr = zeros((nz,nbin,len(wavelength)))
od_mn = zeros((nz,nbin,len(wavelength)))
od_na = zeros((nz,nbin,len(wavelength)))
od_zn = zeros((nz,nbin,len(wavelength)))
od_k = zeros((nz,nbin,len(wavelength)))
od_al = zeros((nz,nbin,len(wavelength)))

od_tot = zeros((nz,nbin,len(wavelength)))

frac_ti = zeros((nz,nbin,len(wavelength)))
frac_mg = zeros((nz,nbin,len(wavelength)))
frac_fe = zeros((nz,nbin,len(wavelength)))
frac_cr = zeros((nz,nbin,len(wavelength)))
frac_mn = zeros((nz,nbin,len(wavelength)))
frac_na = zeros((nz,nbin,len(wavelength)))
frac_zn = zeros((nz,nbin,len(wavelength)))
frac_k = zeros((nz,nbin,len(wavelength)))
frac_al = zeros((nz,nbin,len(wavelength)))

ssa_tot = zeros((nz,nbin,len(wavelength)))

g_tot = zeros((nz,nbin,len(wavelength)))


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

sum_g = zeros((nz,len(wavelength)))
sum_ssa = zeros((nz,len(wavelength)))
	

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
			od_ti[k,j,idw] = t_tio2[k,j]*extinct_ti[idw,j]
			od_mg[k,j,idw] = t_het_mg2sio4[k,j]*extinct_mg[idw,j]
			od_mn[k,j,idw] = t_mns[k,j]*extinct_mn[idw,j]
			od_na[k,j,idw] = t_na2s[k,j]*extinct_na[idw,j]
			od_zn[k,j,idw] = t_zns[k,j]*extinct_zn[idw,j]
			od_k[k,j,idw] = t_kcl[k,j]*extinct_k[idw,j]
			od_al[k,j,idw] = t_al2o3[k,j]*extinct_al[idw,j]
			od_fe[k,j,idw] = (t_het_fe[k,j]+t_pure_fe[k,j])*extinct_fe[idw,j]
			od_cr[k,j,idw] = (t_cr[k,j]+t_het_cr[k,j])*extinct_cr[idw,j]
			
			od_tot[k,j,idw] = od_ti[k,j,idw]+od_mg[k,j,idw]+od_mn[k,j,idw]+od_na[k,j,idw]+od_zn[k,j,idw]+od_k[k,j,idw]+od_al[k,j,idw]+od_fe[k,j,idw]+od_cr[k,j,idw]
			
		
		sum_od_ti[k,idw] = sum(od_ti[k,:,idw])	
		sum_od_mg[k,idw] = sum(od_mg[k,:,idw])	
		sum_od_mn[k,idw] = sum(od_mn[k,:,idw])	
		sum_od_na[k,idw] = sum(od_na[k,:,idw])	
		sum_od_zn[k,idw] = sum(od_zn[k,:,idw])	
		sum_od_k[k,idw] = sum(od_k[k,:,idw])	
		sum_od_al[k,idw] = sum(od_al[k,:,idw])	
		sum_od_fe[k,idw] = sum(od_fe[k,:,idw])	
		sum_od_cr[k,idw] = sum(od_cr[k,:,idw])	
		
		
		
		tot_cloud_opaque[k,idw] = sum_od_ti[k,idw]+sum_od_mg[k,idw]+sum_od_mn[k,idw]+sum_od_na[k,idw]+sum_od_zn[k,idw]+sum_od_k[k,idw]+sum_od_al[k,idw]+sum_od_fe[k,idw]+sum_od_cr[k,idw]


do_Natasha = 0

if do_Natasha:
	print("time for SSA and g analysis..!")
	for idw,w in enumerate(wavelength):
		for k in range(nz): 
			for j in range(nbin):
				if tot_cloud_opaque[k,idw] > 0.:
					frac_ti[k,j,idw] = od_ti[k,j,idw]/tot_cloud_opaque[k,idw]
					frac_mg[k,j,idw] = od_mg[k,j,idw]/tot_cloud_opaque[k,idw]
					frac_mn[k,j,idw] = od_mn[k,j,idw]/tot_cloud_opaque[k,idw]
					frac_na[k,j,idw] = od_na[k,j,idw]/tot_cloud_opaque[k,idw]
					frac_zn[k,j,idw] = od_zn[k,j,idw]/tot_cloud_opaque[k,idw]
					frac_k[k,j,idw] = od_k[k,j,idw]/tot_cloud_opaque[k,idw]
					frac_al[k,j,idw] = od_al[k,j,idw]/tot_cloud_opaque[k,idw]
					frac_fe[k,j,idw] = od_fe[k,j,idw]/tot_cloud_opaque[k,idw]
					frac_cr[k,j,idw] = od_cr[k,j,idw]/tot_cloud_opaque[k,idw]
				'''
				else:
					frac_ti[k,j,idw] = 0.
					frac_mg[k,j,idw] =0.
					frac_mn[k,j,idw] =0.
					frac_na[k,j,idw] =0.
					frac_zn[k,j,idw] =0.
					frac_k[k,j,idw] = 0.
					frac_al[k,j,idw] =0.
					frac_fe[k,j,idw] =0.
					frac_cr[k,j,idw] =0.
				'''

			
				g_tot[k,j,idw] = (frac_ti[k,j,idw]*g_ti[idw,j])+(frac_mg[k,j,idw]*g_mg[idw,j])+(frac_mn[k,j,idw]*g_mn[idw,j])+(frac_na[k,j,idw]*g_na[idw,j])+(
								  frac_zn[k,j,idw]*g_zn[idw,j])+(frac_k[k,j,idw]*g_k[idw,j])+(frac_al[k,j,idw]*g_al[idw,j])+(frac_fe[k,j,idw]*g_fe[idw,j])+(
								  frac_cr[k,j,idw]*g_cr[idw,j])
				ssa_tot[k,j,idw] = (frac_ti[k,j,idw]*ssa_ti[idw,j])+(frac_mg[k,j,idw]*ssa_mg[idw,j])+(frac_mn[k,j,idw]*ssa_mn[idw,j])+(frac_na[k,j,idw]*ssa_na[idw,j])+(
				                    frac_zn[k,j,idw]*ssa_zn[idw,j])+(frac_k[k,j,idw]*ssa_k[idw,j])+(frac_al[k,j,idw]*ssa_al[idw,j])+(frac_fe[k,j,idw]*ssa_fe[idw,j])+(
				                    frac_cr[k,j,idw]*ssa_cr[idw,j])
			sum_g[k,idw] = sum(g_tot[k,:,idw])
			sum_ssa[k,idw] = sum(ssa_tot[k,:,idw])
	
	print('save files')
	np.savetxt(path+'/cloud_ssa_'+name+'.txt', sum_ssa, fmt='%.3e')
	np.savetxt(path+'/cloud_g_'+name+'.txt', sum_g, fmt='%.3e')
	np.savetxt(path+'/cloud_op_'+name+'.txt', tot_cloud_opaque, fmt='%.3e')
	combined = np.vstack((p/1.e6,t,z))
	np.savetxt(path+'/planet_props_'+name+'.txt', combined.T, fmt='%.3e')
	np.savetxt(path+'/wavelength_grid.txt', wavelength.T, fmt='%.3e')


f_tot_opaque = interp2d(wavelength,p/1.e6,tot_cloud_opaque)
i,pressure,temperature = np.loadtxt(tp_file, unpack=True, skiprows = 1)
pressure*=1e-5

print(p/1.e6, pressure)

output = f_tot_opaque(wavelength,pressure)
output *= 100. #opacity is in units of m^-1

np.savetxt(path+'/cloud_opacities_'+name+'.txt', output.T, fmt='%.3e')

plotty = 0

if plotty:
	print("once you plot you don't stop")
	X,Y = np.meshgrid(wavelength, pressure)
	f,ax = plt.subplots(1,1)
	ax.set_axis_bgcolor('k')
	im = ax.pcolormesh(X,Y,output, cmap='Blues', norm=LogNorm())
	im.set_clim(vmin=1.e-19,vmax=1.e1)
	cbar = plt.colorbar(im)
	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.set_xlim(min(wavelength),max(wavelength))
	ax.set_ylim(max(pressure[::-1]),min(pressure[::-1]))
	ax.set_ylabel(r'Pressure [bar]')
	ax.set_xlabel(r'Wavelength [$\mu$m]')
	plt.savefig(path+'/cloud_opacities_'+name+'.jpg', format = 'jpg', bbox_inches='tight')
	plt.close()



		
			