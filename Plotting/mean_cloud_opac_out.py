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

#./mean_cloud_opac_out.py ../../../New_ExoCARMA/apr7_half_latentheat/run/carma/hjt_west_1800.txt ../T_P/t_p_1800_west_limb.dat 1800_west apr11_halflatent_cross_opacities/ > log.txt &

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
print 'Read in output file'
data=open(file,'r')
line = data.next().split()

nz,ngroup,nelem,nbin,ngas,nstep,iskip = map(int,line)
r,ms,dr,rl,ru = zeros((nbin,ngroup)),zeros((nbin,ngroup)),zeros((nbin,ngroup)),zeros((nbin,ngroup)),zeros((nbin,ngroup))
p,t,z = zeros(nz),zeros(nz),zeros(nz)

ntime=nstep/iskip

nz,ngroup,nelem,nbin,ngas,nstep,iskip = map(int,line)
r,ms,dr,rl,ru = zeros((nbin,ngroup)),zeros((nbin,ngroup)),zeros((nbin,ngroup)),zeros((nbin,ngroup)),zeros((nbin,ngroup))
p,t,z = zeros(nz),zeros(nz),zeros(nz)
time = zeros(ntime)

tio2, het_fe, het_mg2sio4, pure_fe, het_cr, cr, mns, na2s,zns,kcl,al2o3 = zeros((nz,nbin,ntime)),zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime))
mass_tio2, mass_het_fe, mass_het_mg2sio4, mass_pure_fe, mass_het_cr, mass_cr, mass_mns, mass_na2s, mass_zns, mass_kcl, mass_al2o3 = zeros((nz,nbin,ntime)),zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime))
t_tio2, t_het_fe, t_het_mg2sio4, t_pure_fe, t_het_cr, t_cr, t_mns, t_na2s,t_zns,t_kcl,t_al2o3 = zeros((nz,nbin)),zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin))
t_mass_tio2, t_mass_het_fe, t_mass_het_mg2sio4, t_mass_pure_fe, t_mass_het_cr, t_mass_cr, t_mass_mns, t_mass_na2s, t_mass_zns, t_mass_kcl, t_mass_al2o3 = zeros((nz,nbin)),zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin))

cross_ti, cross_fe, cross_mg, cross_cr, cross_mn, cross_na, cross_zn, cross_kc, cross_al = zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime))
area_ti, area_fe, area_mg, area_cr, area_mn, area_na, area_zn, area_kc, area_al = zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime))
cross_ti_tot, cross_fe_tot, cross_mg_tot, cross_cr_tot, cross_mn_tot, cross_na_tot, cross_zn_tot, cross_kc_tot, cross_al_tot = zeros((nz,nbin)),zeros((nz,nbin)),zeros((nz,nbin)),zeros((nz,nbin)),zeros((nz,nbin)),zeros((nz,nbin)),zeros((nz,nbin)),zeros((nz,nbin)),zeros((nz,nbin)) 
area_ti_tot, area_fe_tot, area_mg_tot, area_cr_tot, area_mn_tot, area_na_tot, area_zn_tot, area_kc_tot, area_al_tot = zeros((nz,nbin)),zeros((nz,nbin)),zeros((nz,nbin)),zeros((nz,nbin)),zeros((nz,nbin)),zeros((nz,nbin)),zeros((nz,nbin)),zeros((nz,nbin)),zeros((nz,nbin)) 


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
		line = data.next().split()
		r[j,i] = float(line[2])
		ms[j,i] = float(line[3])
		dr[j,i] = float(line[4])
		rl[j,i] = float(line[5])
		ru[j,i] = float(line[6])
for i in range(nz):
	line = data.next().split()
	z[i] = float(line[1])
	p[i] = float(line[3])
	t[i] = float(line[4])
			
for i in range(ntime-1):
	time[i] = float(data.next())
	for j in range(nbin):
		for k in range(nz):
			line = data.next().split()
			
			tio2[k,j,i] = float(line[2])
			mass_tio2[k,j,i] = tio2[k,j,i]*(M_TIO2*((4./3.)*np.pi*(r[j,0]*1e-4)**3.))
			
			al2o3[k,j,i] = float(line[3])
			mass_al2o3[k,j,i] = al2o3[k,j,i]*(M_AL2O3*((4./3.)*np.pi*(r[j,1]*1e-4)**3.))
			
			het_fe[k,j,i] = float(line[4])
			mass_het_fe[k,j,i] = het_fe[k,j,i]*(M_FE*((4./3.)*np.pi*(r[j,2]*1e-4)**3.))
			
			het_mg2sio4[k,j,i] = float(line[5])
			mass_het_mg2sio4[k,j,i] = het_mg2sio4[k,j,i]*(M_MG2SIO4*((4./3.)*np.pi*(r[j,3]*1e-4)**3.))
			
			het_cr[k,j,i] = float(line[6])
			mass_het_cr[k,j,i] = het_cr[k,j,i]*(M_CR*((4./3.)*np.pi*(r[j,4]*1e-4)**3.))
			
			mns[k,j,i] = float(line[7])#/(math.log(ru[j,5]) - math.log(rl[j,5]))
			mass_mns[k,j,i] = mns[k,j,i]*(M_MNS*((4./3.)*np.pi*(r[j,5]*1e-4)**3.))
			
			na2s[k,j,i] = float(line[8])#/(math.log(ru[j,6]) - math.log(rl[j,6]))
			mass_na2s[k,j,i] = na2s[k,j,i]*(M_NA2S*((4./3.)*np.pi*(r[j,6]*1e-4)**3.))
			
			pure_fe[k,j,i] = float(line[9])#/(math.log(ru[j,7]) - math.log(rl[j,7]))
			mass_pure_fe[k,j,i] = pure_fe[k,j,i]*(M_FE*((4./3.)*np.pi*(r[j,7]*1e-4)**3.))
			
			cr[k,j,i] = float(line[10])#/(math.log(ru[j,8]) - math.log(rl[j,8]))
			mass_cr[k,j,i] = cr[k,j,i]*(M_CR*((4./3.)*np.pi*(r[j,8]*1e-4)**3.))
			
			kcl[k,j,i] = float(line[11])#/(math.log(ru[j,9]) - math.log(rl[j,9]))
			mass_kcl[k,j,i] = kcl[k,j,i]*(M_KCL*((4./3.)*np.pi*(r[j,9]*1e-4)**3.))
			
			zns[k,j,i] = float(line[12])#/(math.log(ru[j,10]) - math.log(rl[j,10]))
			mass_zns[k,j,i] = zns[k,j,i]*(M_ZNS*((4./3.)*np.pi*(r[j,10]*1e-4)**3.))
			
			cross_ti[k,j,i] = tio2[k,j,i]*np.pi*(r[j,0]*1e-4)**2.
			cross_fe[k,j,i] = (het_fe[k,j,i]*np.pi*(r[j,2]*1e-4)**2.)+(pure_fe[k,j,i]*np.pi*(r[j,7]*1e-4)**2.)
			cross_mg[k,j,i] = (het_mg2sio4[k,j,i])*np.pi*(r[j,3]*1e-4)**2.
			cross_cr[k,j,i] = (het_cr[k,j,i]*np.pi*(r[j,4]*1e-4)**2.)+(cr[k,j,i]*np.pi*(r[j,8]*1e-4)**2.)
			cross_mn[k,j,i] = mns[k,j,i]*np.pi*(r[j,5]*1e-4)**2.
			cross_na[k,j,i] = na2s[k,j,i]*np.pi*(r[j,6]*1e-4)**2.
			cross_zn[k,j,i] = zns[k,j,i]*np.pi*(r[j,10]*1e-4)**2.
			cross_kc[k,j,i] = kcl[k,j,i]*np.pi*(r[j,9]*1e-4)**2.
			cross_al[k,j,i] = al2o3[k,j,i]*np.pi*(r[j,1]*1e-4)**2.
			
			area_ti[k,j,i] = tio2[k,j,i]*np.pi*(r[j,0]*1e-4)**3.
			area_fe[k,j,i] = (het_fe[k,j,i]*np.pi*(r[j,2]*1e-4)**3.)+(pure_fe[k,j,i]*np.pi*(r[j,7]*1e-4)**3.)
			area_mg[k,j,i] = (het_mg2sio4[k,j,i])*np.pi*(r[j,3]*1e-4)**3.
			area_cr[k,j,i] = (het_cr[k,j,i]*np.pi*(r[j,4]*1e-4)**3.)+(cr[k,j,i]*np.pi*(r[j,8]*1e-4)**3.)
			area_mn[k,j,i] = mns[k,j,i]*np.pi*(r[j,5]*1e-4)**3.
			area_na[k,j,i] = na2s[k,j,i]*np.pi*(r[j,6]*1e-4)**3.
			area_zn[k,j,i] = zns[k,j,i]*np.pi*(r[j,10]*1e-4)**3.
			area_kc[k,j,i] = kcl[k,j,i]*np.pi*(r[j,9]*1e-4)**3.
			area_al[k,j,i] = al2o3[k,j,i]*np.pi*(r[j,1]*1e-4)**3.
			
	else:
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
			
			cross_ti_tot[:,:] += cross_ti[:,:,i]
			cross_mg_tot[:,:] += cross_mg[:,:,i]
			cross_fe_tot[:,:] += cross_fe[:,:,i]
			cross_cr_tot[:,:] += cross_cr[:,:,i]
			cross_zn_tot[:,:] += cross_zn[:,:,i]
			cross_na_tot[:,:] += cross_na[:,:,i]
			cross_al_tot[:,:] += cross_al[:,:,i]
			cross_kc_tot[:,:] += cross_kc[:,:,i]
			cross_mn_tot[:,:] += cross_mn[:,:,i]
			
			area_ti_tot[:,:] += area_ti[:,:,i]
			area_mg_tot[:,:] += area_mg[:,:,i]
			area_fe_tot[:,:] += area_fe[:,:,i]
			area_cr_tot[:,:] += area_cr[:,:,i]
			area_zn_tot[:,:] += area_zn[:,:,i]
			area_na_tot[:,:] += area_na[:,:,i]
			area_al_tot[:,:] += area_al[:,:,i]
			area_kc_tot[:,:] += area_kc[:,:,i]
			area_mn_tot[:,:] += area_mn[:,:,i]
			
						
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

cross_ti_tot /= count
cross_mg_tot /= count
cross_fe_tot /= count
cross_cr_tot /= count
cross_al_tot /= count
cross_mn_tot /= count
cross_na_tot /= count
cross_zn_tot /= count
cross_kc_tot /= count

area_ti_tot /= count
area_mg_tot /= count
area_fe_tot /= count
area_cr_tot /= count
area_al_tot /= count
area_mn_tot /= count
area_na_tot /= count
area_zn_tot /= count
area_kc_tot /= count



num_density_ti = zeros((nz))
num_density_mg = zeros((nz))
num_density_fe = zeros((nz))
num_density_cr = zeros((nz))
num_density_mn = zeros((nz))
num_density_na = zeros((nz))
num_density_zn = zeros((nz))
num_density_kc = zeros((nz))
num_density_al = zeros((nz))

ps_ti_mass = zeros((nz))
ps_mg_mass = zeros((nz))
ps_fe_mass = zeros((nz))
ps_cr_mass = zeros((nz))
ps_mn_mass = zeros((nz))
ps_na_mass = zeros((nz))
ps_zn_mass = zeros((nz))
ps_kc_mass = zeros((nz))
ps_al_mass = zeros((nz))

ps_ti_area = zeros((nz))
ps_mg_area = zeros((nz))
ps_fe_area = zeros((nz))
ps_cr_area = zeros((nz))
ps_mn_area = zeros((nz))
ps_na_area = zeros((nz))
ps_zn_area = zeros((nz))
ps_kc_area = zeros((nz))
ps_al_area = zeros((nz))

ps_ti_cross = zeros((nz))
ps_mg_cross = zeros((nz))
ps_fe_cross = zeros((nz))
ps_cr_cross = zeros((nz))
ps_mn_cross = zeros((nz))
ps_na_cross = zeros((nz))
ps_zn_cross = zeros((nz))
ps_kc_cross = zeros((nz))
ps_al_cross = zeros((nz))


numdens_ti_mass = zeros((nz))
numdens_mg_mass = zeros((nz))
numdens_fe_mass = zeros((nz))
numdens_cr_mass = zeros((nz))
numdens_mn_mass = zeros((nz))
numdens_na_mass = zeros((nz))
numdens_zn_mass = zeros((nz))
numdens_kc_mass = zeros((nz))
numdens_al_mass = zeros((nz))

numdens_ti_cross = zeros((nz))
numdens_mg_cross = zeros((nz))
numdens_fe_cross = zeros((nz))
numdens_cr_cross = zeros((nz))
numdens_mn_cross = zeros((nz))
numdens_na_cross = zeros((nz))
numdens_zn_cross = zeros((nz))
numdens_kc_cross = zeros((nz))
numdens_al_cross = zeros((nz))

numdens_ti_area = zeros((nz))
numdens_mg_area = zeros((nz))
numdens_fe_area = zeros((nz))
numdens_cr_area = zeros((nz))
numdens_mn_area = zeros((nz))
numdens_na_area = zeros((nz))
numdens_zn_area = zeros((nz))
numdens_kc_area = zeros((nz))
numdens_al_area = zeros((nz))

for k in range(nz):
	num_density_ti[k] = np.sum(t_tio2[k,:])
	num_density_mg[k] = np.sum(t_het_mg2sio4[k,:])
	num_density_fe[k] = np.sum(t_het_fe[k,:]+t_pure_fe[k,:])
	num_density_cr[k] = np.sum(t_cr[k,:]+t_het_cr[k,:])
	num_density_mn[k] = np.sum(t_mns[k,:])
	num_density_na[k] = np.sum(t_na2s[k,:])
	num_density_zn[k] = np.sum(t_zns[k,:])
	num_density_kc[k] = np.sum(t_kcl[k,:])
	num_density_al[k] = np.sum(t_al2o3[k,:])
		
	ps_ti_mass[k] = np.sum(t_mass_tio2[k,:])/(num_density_ti[k])
	ps_ti_mass[k] = (ps_ti_mass[k]/M_TIO2/(4./3.)/np.pi)**(1./3.)*1.e4 #microns	
	ps_mg_mass[k] = np.sum(t_mass_het_mg2sio4[k,:])/(num_density_mg[k])
	ps_mg_mass[k] = (ps_mg_mass[k]/M_MG2SIO4/(4./3.)/np.pi)**(1./3.)*1.e4
	ps_fe_mass[k] = np.sum(t_mass_het_fe[k,:]+t_mass_pure_fe[k,:])/(num_density_fe[k])
	ps_fe_mass[k] = (ps_fe_mass[k]/M_FE/(4./3.)/np.pi)**(1./3.)*1.e4
	ps_cr_mass[k] = np.sum(t_mass_cr[k,:]+t_mass_het_cr[k,:])/(num_density_cr[k])
	ps_cr_mass[k] = (ps_cr_mass[k]/M_CR/(4./3.)/np.pi)**(1./3.)*1.e4
	ps_mn_mass[k] = np.sum(t_mass_mns[k,:])/(num_density_mn[k])
	ps_mn_mass[k] = (ps_mn_mass[k]/M_MNS/(4./3.)/np.pi)**(1./3.)*1.e4
	ps_na_mass[k] = np.sum(t_mass_na2s[k,:])/(num_density_na[k])
	ps_na_mass[k] = (ps_na_mass[k]/M_NA2S/(4./3.)/np.pi)**(1./3.)*1.e4
	ps_zn_mass[k] = np.sum(t_mass_zns[k,:])/(num_density_zn[k])
	ps_zn_mass[k] = (ps_zn_mass[k]/M_ZNS/(4./3.)/np.pi)**(1./3.)*1.e4
	ps_kc_mass[k] = np.sum(t_mass_kcl[k,:])/(num_density_kc[k])
	ps_kc_mass[k] = (ps_kc_mass[k]/M_KCL/(4./3.)/np.pi)**(1./3.)*1.e4
	ps_al_mass[k] = np.sum(t_mass_al2o3[k,:])/(num_density_al[k])
	ps_al_mass[k] = (ps_al_mass[k]/M_AL2O3/(4./3.)/np.pi)**(1./3.)*1.e4
	
	numdens_ti_mass[k] = np.sum(t_mass_tio2[k,:])/(M_TIO2*(4./3.)*np.pi*(ps_ti_mass[k]*1.e-4)**3.)
	numdens_mg_mass[k] = np.sum(t_mass_het_mg2sio4[k,:])/(M_MG2SIO4*(4./3.)*np.pi*(ps_mg_mass[k]*1.e-4)**3.)
	numdens_fe_mass[k] = np.sum(t_mass_het_fe[k,:]+t_mass_pure_fe[k,:])/(M_FE*(4./3.)*np.pi*(ps_fe_mass[k]*1.e-4)**3.)
	numdens_cr_mass[k] = np.sum(t_mass_cr[k,:]+t_mass_het_cr[k,:])/(M_CR*(4./3.)*np.pi*(ps_cr_mass[k]*1.e-4)**3.)
	numdens_mn_mass[k] = np.sum(t_mass_mns[k,:])/(M_MNS*(4./3.)*np.pi*(ps_mn_mass[k]*1.e-4)**3.)
	numdens_na_mass[k] = np.sum(t_mass_na2s[k,:])/(M_NA2S*(4./3.)*np.pi*(ps_na_mass[k]*1.e-4)**3.)
	numdens_zn_mass[k] = np.sum(t_mass_zns[k,:])/(M_ZNS*(4./3.)*np.pi*(ps_zn_mass[k]*1.e-4)**3.)
	numdens_kc_mass[k] = np.sum(t_mass_kcl[k,:])/(M_KCL*(4./3.)*np.pi*(ps_kc_mass[k]*1.e-4)**3.)
	numdens_al_mass[k] = np.sum(t_mass_al2o3[k,:])/(M_AL2O3*(4./3.)*np.pi*(ps_al_mass[k]*1.e-4)**3.)
	
	ps_ti_cross[k] = np.sum(cross_ti_tot[k,:])/(num_density_ti[k])
	ps_ti_cross[k] = (ps_ti_cross[k]/np.pi)**.5*1.e4 #microns	
	ps_mg_cross[k] = np.sum(cross_mg_tot[k,:])/(num_density_mg[k])
	ps_mg_cross[k] = (ps_mg_cross[k]/np.pi)**.5*1.e4
	ps_fe_cross[k] = np.sum(cross_fe_tot[k,:])/(num_density_fe[k])
	ps_fe_cross[k] = (ps_fe_cross[k]/np.pi)**.5*1.e4
	ps_cr_cross[k] = np.sum(cross_cr_tot[k,:])/(num_density_cr[k])
	ps_cr_cross[k] = (ps_cr_cross[k]/np.pi)**.5*1.e4
	ps_mn_cross[k] = np.sum(cross_mn_tot[k,:])/(num_density_mn[k])
	ps_mn_cross[k] = (ps_mn_cross[k]/np.pi)**.5*1.e4
	ps_na_cross[k] = np.sum(cross_na_tot[k,:])/(num_density_na[k])
	ps_na_cross[k] = (ps_na_cross[k]/np.pi)**.5*1.e4
	ps_zn_cross[k] = np.sum(cross_zn_tot[k,:])/(num_density_zn[k])
	ps_zn_cross[k] = (ps_zn_cross[k]/np.pi)**.5*1.e4
	ps_kc_cross[k] = np.sum(cross_kc_tot[k,:])/(num_density_kc[k])
	ps_kc_cross[k] = (ps_kc_cross[k]/np.pi)**.5*1.e4
	ps_al_cross[k] = np.sum(cross_al_tot[k,:])/(num_density_al[k])
	ps_al_cross[k] = (ps_al_cross[k]/np.pi)**.5*1.e4
	
	numdens_ti_cross[k] = np.sum(t_mass_tio2[k,:])/(M_TIO2*(4./3.)*np.pi*(ps_ti_cross[k]*1.e-4)**3.)
	numdens_mg_cross[k] = np.sum(t_mass_het_mg2sio4[k,:])/(M_MG2SIO4*(4./3.)*np.pi*(ps_mg_cross[k]*1.e-4)**3.)
	numdens_fe_cross[k] = np.sum(t_mass_het_fe[k,:]+t_mass_pure_fe[k,:])/(M_FE*(4./3.)*np.pi*(ps_fe_cross[k]*1.e-4)**3.)
	numdens_cr_cross[k] = np.sum(t_mass_cr[k,:]+t_mass_het_cr[k,:])/(M_CR*(4./3.)*np.pi*(ps_cr_cross[k]*1.e-4)**3.)
	numdens_mn_cross[k] = np.sum(t_mass_mns[k,:])/(M_MNS*(4./3.)*np.pi*(ps_mn_cross[k]*1.e-4)**3.)
	numdens_na_cross[k] = np.sum(t_mass_na2s[k,:])/(M_NA2S*(4./3.)*np.pi*(ps_na_cross[k]*1.e-4)**3.)
	numdens_zn_cross[k] = np.sum(t_mass_zns[k,:])/(M_ZNS*(4./3.)*np.pi*(ps_zn_cross[k]*1.e-4)**3.)
	numdens_kc_cross[k] = np.sum(t_mass_kcl[k,:])/(M_KCL*(4./3.)*np.pi*(ps_kc_cross[k]*1.e-4)**3.)
	numdens_al_cross[k] = np.sum(t_mass_al2o3[k,:])/(M_AL2O3*(4./3.)*np.pi*(ps_al_cross[k]*1.e-4)**3.)
	
	ps_ti_area[k] = np.sum(area_ti_tot[k,:])/np.sum(cross_ti_tot[k,:])
	ps_ti_area[k] = ps_ti_area[k]*1.e4 #microns	
	ps_mg_area[k] = np.sum(area_mg_tot[k,:])/np.sum(cross_mg_tot[k,:])
	ps_mg_area[k] = ps_mg_area[k]*1.e4
	ps_fe_area[k] = np.sum(area_fe_tot[k,:])/np.sum(cross_fe_tot[k,:])
	ps_fe_area[k] = ps_fe_area[k]*1.e4
	ps_cr_area[k] = np.sum(area_cr_tot[k,:])/np.sum(cross_cr_tot[k,:])
	ps_cr_area[k] = ps_cr_area[k]*1.e4
	ps_mn_area[k] = np.sum(area_mn_tot[k,:])/np.sum(cross_mn_tot[k,:])
	ps_mn_area[k] = ps_mn_area[k]*1.e4
	ps_na_area[k] = np.sum(area_na_tot[k,:])/np.sum(cross_na_tot[k,:])
	ps_na_area[k] = ps_na_area[k]*1.e4
	ps_zn_area[k] = np.sum(area_zn_tot[k,:])/np.sum(cross_zn_tot[k,:])
	ps_zn_area[k] = ps_zn_area[k]*1.e4
	ps_kc_area[k] = np.sum(area_kc_tot[k,:])/np.sum(cross_kc_tot[k,:])
	ps_kc_area[k] = ps_kc_area[k]*1.e4
	ps_al_area[k] = np.sum(area_al_tot[k,:])/np.sum(cross_al_tot[k,:])
	ps_al_area[k] = ps_al_area[k]*1.e4
	
	numdens_ti_area[k] = np.sum(t_mass_tio2[k,:])/(M_TIO2*(4./3.)*np.pi*(ps_ti_area[k]*1.e-4)**3.)
	numdens_mg_area[k] = np.sum(t_mass_het_mg2sio4[k,:])/(M_MG2SIO4*(4./3.)*np.pi*(ps_mg_area[k]*1.e-4)**3.)
	numdens_fe_area[k] = np.sum(t_mass_het_fe[k,:]+t_mass_pure_fe[k,:])/(M_FE*(4./3.)*np.pi*(ps_fe_area[k]*1.e-4)**3.)
	numdens_cr_area[k] = np.sum(t_mass_cr[k,:]+t_mass_het_cr[k,:])/(M_CR*(4./3.)*np.pi*(ps_cr_area[k]*1.e-4)**3.)
	numdens_mn_area[k] = np.sum(t_mass_mns[k,:])/(M_MNS*(4./3.)*np.pi*(ps_mn_area[k]*1.e-4)**3.)
	numdens_na_area[k] = np.sum(t_mass_na2s[k,:])/(M_NA2S*(4./3.)*np.pi*(ps_na_area[k]*1.e-4)**3.)
	numdens_zn_area[k] = np.sum(t_mass_zns[k,:])/(M_ZNS*(4./3.)*np.pi*(ps_zn_area[k]*1.e-4)**3.)
	numdens_kc_area[k] = np.sum(t_mass_kcl[k,:])/(M_KCL*(4./3.)*np.pi*(ps_kc_area[k]*1.e-4)**3.)
	numdens_al_area[k] = np.sum(t_mass_al2o3[k,:])/(M_AL2O3*(4./3.)*np.pi*(ps_al_area[k]*1.e-4)**3.)


############
# Read in mie scattering output
############

wavelength,extinct_ti = np.loadtxt('rf_data/ti_rad_1.txt', unpack=True)


od_ti_m = zeros((nz,len(wavelength)))
od_mg_m = zeros((nz,len(wavelength)))
od_fe_m = zeros((nz,len(wavelength)))
od_cr_m = zeros((nz,len(wavelength)))
od_mn_m = zeros((nz,len(wavelength)))
od_na_m = zeros((nz,len(wavelength)))
od_zn_m = zeros((nz,len(wavelength)))
od_k_m = zeros((nz,len(wavelength)))
od_al_m = zeros((nz,len(wavelength)))
tot_cloud_opaque_m = zeros((nz,len(wavelength)))

od_ti_a = zeros((nz,len(wavelength)))
od_mg_a = zeros((nz,len(wavelength)))
od_fe_a = zeros((nz,len(wavelength)))
od_cr_a = zeros((nz,len(wavelength)))
od_mn_a = zeros((nz,len(wavelength)))
od_na_a = zeros((nz,len(wavelength)))
od_zn_a = zeros((nz,len(wavelength)))
od_k_a = zeros((nz,len(wavelength)))
od_al_a = zeros((nz,len(wavelength)))
tot_cloud_opaque_a = zeros((nz,len(wavelength)))

od_ti_c = zeros((nz,len(wavelength)))
od_mg_c = zeros((nz,len(wavelength)))
od_fe_c = zeros((nz,len(wavelength)))
od_cr_c = zeros((nz,len(wavelength)))
od_mn_c = zeros((nz,len(wavelength)))
od_na_c = zeros((nz,len(wavelength)))
od_zn_c = zeros((nz,len(wavelength)))
od_k_c = zeros((nz,len(wavelength)))
od_al_c = zeros((nz,len(wavelength)))
tot_cloud_opaque_c = zeros((nz,len(wavelength)))
	

extinct_ti = zeros((len(wavelength),nbin))
extinct_mg = zeros((len(wavelength),nbin))
extinct_fe = zeros((len(wavelength),nbin))
extinct_cr = zeros((len(wavelength),nbin))
extinct_mn = zeros((len(wavelength),nbin))
extinct_na = zeros((len(wavelength),nbin))
extinct_zn = zeros((len(wavelength),nbin))
extinct_k = zeros((len(wavelength),nbin))
extinct_al = zeros((len(wavelength),nbin))

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
		extinct_ti[idw,j] = ext_tio2[count]
		extinct_mg[idw,j] = ext_mg2sio4[count]
		extinct_mn[idw,j] = ext_mns[count]
		extinct_na[idw,j] = ext_na2s[count]
		extinct_k[idw,j] = ext_kcl[count]
		extinct_fe[idw,j] = ext_fe[count]
		extinct_cr[idw,j] = ext_cr[count]
		extinct_al[idw,j] = ext_al2o3[count]
		extinct_zn[idw,j] = ext_zns[count]
		count +=1

f_ex_ti = interp2d(r[:,1],wavelength,extinct_ti, fill_value = 0.)
f_ex_mg = interp2d(r[:,1],wavelength,extinct_mg, fill_value = 0.)
f_ex_mn = interp2d(r[:,1],wavelength,extinct_mn, fill_value = 0.)
f_ex_na = interp2d(r[:,1],wavelength,extinct_na, fill_value = 0.)
f_ex_zn = interp2d(r[:,1],wavelength,extinct_zn, fill_value = 0.)
f_ex_al = interp2d(r[:,1],wavelength,extinct_al, fill_value = 0.)
f_ex_kc = interp2d(r[:,1],wavelength,extinct_k, fill_value = 0.)
f_ex_cr = interp2d(r[:,1],wavelength,extinct_cr, fill_value = 0.)
f_ex_fe = interp2d(r[:,1],wavelength,extinct_fe, fill_value = 0.)
	
##########################################################################################
##########################################################################################	
print 'BIG LOOP TIME!'
##########################################################################################
########################################################################################## 

for idw,w in enumerate(wavelength):
	for k in range(nz):
		if np.isnan(ps_ti_mass[k]) == False: 
			od_ti_m[k,idw] = f_ex_ti(ps_ti_mass[k],w)[0]*numdens_ti_mass[k]
		if np.isnan(ps_mg_mass[k]) == False: 
			od_mg_m[k,idw] = f_ex_ti(ps_mg_mass[k],w)[0]*numdens_mg_mass[k]
		if np.isnan(ps_mn_mass[k]) == False: 
			od_mn_m[k,idw] = f_ex_ti(ps_mn_mass[k],w)[0]*numdens_mn_mass[k] 
		if np.isnan(ps_na_mass[k]) == False: 
			od_na_m[k,idw] = f_ex_ti(ps_na_mass[k],w)[0]*numdens_na_mass[k] 
		if np.isnan(ps_zn_mass[k]) == False: 
			od_zn_m[k,idw] = f_ex_ti(ps_zn_mass[k],w)[0]*numdens_zn_mass[k] 
		if np.isnan(ps_fe_mass[k]) == False: 
			od_fe_m[k,idw] = f_ex_ti(ps_fe_mass[k],w)[0]*numdens_fe_mass[k] 
		if np.isnan(ps_cr_mass[k]) == False: 
			od_cr_m[k,idw] = f_ex_ti(ps_cr_mass[k],w)[0]*numdens_cr_mass[k] 
		if np.isnan(ps_al_mass[k]) == False: 
			od_al_m[k,idw] = f_ex_ti(ps_al_mass[k],w)[0]*numdens_al_mass[k] 
		if np.isnan(ps_kc_mass[k]) == False: 
			od_kc_m[k,idw] = f_ex_ti(ps_kc_mass[k],w)[0]*numdens_kc_mass[k] 
		tot_cloud_opaque_m[k,idw] = od_ti_m[k,idw]+od_mg_m[k,idw]+od_mn_m[k,idw]+od_na_m[k,idw]+od_zn_m[k,idw]+od_k_m[k,idw]+od_al_m[k,idw]+od_fe_m[k,idw]+od_cr_m[k,idw]
		
		if np.isnan(ps_ti_cross[k]) == False: 
			od_ti_c[k,idw] = f_ex_ti(ps_ti_cross[k],w)[0]*numdens_ti_cross[k] 
		if np.isnan(ps_mg_cross[k]) == False: 
			od_mg_c[k,idw] = f_ex_ti(ps_mg_cross[k],w)[0]*numdens_mg_cross[k] 
		if np.isnan(ps_mn_cross[k]) == False: 
			od_mn_c[k,idw] = f_ex_ti(ps_mn_cross[k],w)[0]*numdens_mn_cross[k] 
		if np.isnan(ps_na_cross[k]) == False: 
			od_na_c[k,idw] = f_ex_ti(ps_na_cross[k],w)[0]*numdens_na_cross[k] 
		if np.isnan(ps_zn_cross[k]) == False: 
			od_zn_c[k,idw] = f_ex_ti(ps_zn_cross[k],w)[0]*numdens_zn_cross[k] 
		if np.isnan(ps_fe_cross[k]) == False: 
			od_fe_c[k,idw] = f_ex_ti(ps_fe_cross[k],w)[0]*numdens_fe_cross[k] 
		if np.isnan(ps_cr_cross[k]) == False: 
			od_cr_c[k,idw] = f_ex_ti(ps_cr_cross[k],w)[0]*numdens_cr_cross[k] 
		if np.isnan(ps_al_cross[k]) == False: 
			od_al_c[k,idw] = f_ex_ti(ps_al_cross[k],w)[0]*numdens_al_cross[k] 
		if np.isnan(ps_kc_cross[k]) == False: 
			od_kc_c[k,idw] = f_ex_ti(ps_kc_cross[k],w)[0]*numdens_kc_cross[k] 
		tot_cloud_opaque_c[k,idw] = od_ti_c[k,idw]+od_mg_c[k,idw]+od_mn_c[k,idw]+od_na_c[k,idw]+od_zn_c[k,idw]+od_k_c[k,idw]+od_al_c[k,idw]+od_fe_c[k,idw]+od_cr_c[k,idw]
		
		if np.isnan(ps_ti_area[k]) == False: 
			od_ti_a[k,idw] = f_ex_ti(ps_ti_area[k],w)[0]*numdens_ti_area[k] 
		if np.isnan(ps_mg_area[k]) == False: 
			od_mg_a[k,idw] = f_ex_ti(ps_mg_area[k],w)[0]*numdens_mg_area[k] 
		if np.isnan(ps_mn_area[k]) == False: 
			od_mn_a[k,idw] = f_ex_ti(ps_mn_area[k],w)[0]*numdens_mn_area[k] 
		if np.isnan(ps_na_area[k]) == False: 
			od_na_a[k,idw] = f_ex_ti(ps_na_area[k],w)[0]*numdens_na_area[k] 
		if np.isnan(ps_zn_area[k]) == False: 
			od_zn_a[k,idw] = f_ex_ti(ps_zn_area[k],w)[0]*numdens_zn_area[k] 
		if np.isnan(ps_fe_area[k]) == False: 
			od_fe_a[k,idw] = f_ex_ti(ps_fe_area[k],w)[0]*numdens_fe_area[k] 
		if np.isnan(ps_cr_area[k]) == False: 
			od_cr_a[k,idw] = f_ex_ti(ps_cr_area[k],w)[0]*numdens_cr_area[k] 
		if np.isnan(ps_al_area[k]) == False: 
			od_al_a[k,idw] = f_ex_ti(ps_al_area[k],w)[0]*numdens_al_area[k] 
		if np.isnan(ps_kc_area[k]) == False: 
			od_kc_a[k,idw] = f_ex_ti(ps_kc_area[k],w)[0]*numdens_kc_area[k] 
		tot_cloud_opaque_a[k,idw] = od_ti_a[k,idw]+od_mg_a[k,idw]+od_mn_a[k,idw]+od_na_a[k,idw]+od_zn_a[k,idw]+od_k_a[k,idw]+od_al_a[k,idw]+od_fe_a[k,idw]+od_cr_a[k,idw]


f_tot_opaque_m = interp2d(wavelength,p/1.e6,tot_cloud_opaque_m)
f_tot_opaque_a = interp2d(wavelength,p/1.e6,tot_cloud_opaque_a)
f_tot_opaque_c = interp2d(wavelength,p/1.e6,tot_cloud_opaque_c)
i,pressure,temperature = np.loadtxt(tp_file, unpack=True, skiprows = 1)
pressure*=1e-5

output_m = f_tot_opaque_m(wavelength,pressure[::-1])*100.
np.savetxt(path+'/cloud_opacities_mass_'+name+'.txt', output_m.T, fmt='%.3e')
output_a = f_tot_opaque_a(wavelength,pressure[::-1])*100.
np.savetxt(path+'/cloud_opacities_area_'+name+'.txt', output_a.T, fmt='%.3e')
output_c = f_tot_opaque_c(wavelength,pressure[::-1])*100.
np.savetxt(path+'/cloud_opacities_cross_'+name+'.txt', output_c.T, fmt='%.3e')

plotty = 1

if plotty:
	X,Y = np.meshgrid(wavelength, pressure[::-1])
	f,ax = plt.subplots(1,1)
	ax.set_axis_bgcolor('k')
	im = ax.pcolormesh(X,Y,output_m, cmap='Blues', norm=LogNorm())
	im.set_clim(vmin=1.e-19,vmax=1.e1)
	cbar = plt.colorbar(im)
	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.set_xlim(min(wavelength),max(wavelength))
	ax.set_ylabel(r'Pressure [pascale]')
	ax.set_xlabel(r'Wavelength [$\mu$m]')
	plt.savefig(path+'/co_mass_'+name+'.jpg', format = 'jpg', bbox_inches='tight')
	plt.close()

	X,Y = np.meshgrid(wavelength, pressure[::-1])
	f,ax = plt.subplots(1,1)
	ax.set_axis_bgcolor('k')
	im = ax.pcolormesh(X,Y,output_a, cmap='Blues', norm=LogNorm())
	im.set_clim(vmin=1.e-19,vmax=1.e1)
	cbar = plt.colorbar(im)
	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.set_xlim(min(wavelength),max(wavelength))
	ax.set_ylabel(r'Pressure [pascale]')
	ax.set_xlabel(r'Wavelength [$\mu$m]')
	plt.savefig(path+'/co_area_'+name+'.jpg', format = 'jpg', bbox_inches='tight')
	plt.close()
	
	X,Y = np.meshgrid(wavelength, pressure[::-1])
	f,ax = plt.subplots(1,1)
	ax.set_axis_bgcolor('k')
	im = ax.pcolormesh(X,Y,output_c, cmap='Blues', norm=LogNorm())
	im.set_clim(vmin=1.e-19,vmax=1.e1)
	cbar = plt.colorbar(im)
	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.set_xlim(min(wavelength),max(wavelength))
	ax.set_ylabel(r'Pressure [pascale]')
	ax.set_xlabel(r'Wavelength [$\mu$m]')
	plt.savefig(path+'/co_cross_'+name+'.jpg', format = 'jpg', bbox_inches='tight')
	plt.close()




		
			