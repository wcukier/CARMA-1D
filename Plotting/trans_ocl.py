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

#./cloud_opacity_out_morebins.py hjt_west_1800.txt ../T_P/t_p_1800_west_limb.dat 1800_west opacities/ > log.txt &

#./cloud_opacity_out_morebins.py ../../../New_ExoCARMA/apr7_half_latentheat/run/carma/hjt_west_1800.txt ../T_P/t_p_1800_west_limb.dat 1800_west apr11_halflatent_opacities/ > log.txt &
#./cloud_opacity_out_morebins.py ../../../New_ExoCARMA/apr7_half_latentheat/run/carma/hjt_east_1800.txt ../T_P/t_p_1800_east_limb.dat 1800_east apr11_halflatent_opacities/ > log.txt &
#./cloud_opacity_out_morebins.py ../../../New_ExoCARMA/apr7_half_latentheat/run/carma/hjt_pole_1800.txt ../T_P/t_p_1800_pole.dat 1800_pole apr11_halflatent_opacities/ > log.txt &

# Read in inputs from command line
opts, args = getopt.getopt(sys.argv[1:], "h")
file = args[0]
name = args[1]
path = args[2]

name = str(name)

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


ntime = nstep/iskip

time = zeros(ntime)

tio2, het_fe, het_mg2sio4, pure_fe, het_cr, cr, mns, na2s,zns,kcl,al2o3 = zeros((nz,nbin,ntime)),zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime)), zeros((nz,nbin,ntime))
t_tio2, t_het_fe, t_het_mg2sio4, t_pure_fe, t_het_cr, t_cr, t_mns, t_na2s,t_zns,t_kcl,t_al2o3 = zeros((nz,nbin)),zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin))

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
			
			al2o3[k,j,i] = float(line[3])
			
			het_fe[k,j,i] = float(line[4])
			
			het_mg2sio4[k,j,i] = float(line[5])
			
			het_cr[k,j,i] = float(line[6])
			
			mns[k,j,i] = float(line[7])
			
			na2s[k,j,i] = float(line[8])
			
			pure_fe[k,j,i] = float(line[9])
			
			cr[k,j,i] = float(line[10])
			
			kcl[k,j,i] = float(line[11])
			
			zns[k,j,i] = float(line[12])
	
	if i > ntime-10 :
		count += 1
			
		t_tio2[:,:] += tio2[:,:,i]
			
		t_het_fe[:,:] += het_fe[:,:,i]
			
		t_het_mg2sio4[:,:] += het_mg2sio4[:,:,i]
			
		t_pure_fe[:,:] += pure_fe[:,:,i]
			
		t_het_cr[:,:] += het_cr[:,:,i]
			
		t_cr[:,:] += cr[:,:,i]
			
		t_mns[:,:] += mns[:,:,i]
			
		t_na2s[:,:] += na2s[:,:,i]
			
		t_zns[:,:] += zns[:,:,i]
			
		t_kcl[:,:] += kcl[:,:,i]
			
		t_al2o3[:,:] += al2o3[:,:,i]
		
t_tio2 /= count 
t_het_fe /= count
t_het_mg2sio4 /= count
t_pure_fe /= count
t_het_cr /= count
t_cr /= count
t_mns /= count
t_na2s /= count
t_zns /= count
t_kcl /= count
t_al2o3 /= count

wavelength,extinct_ti = np.loadtxt('/auto/home/dkpowell/EXO_TRANSMIT/UPDATED_Exotransmit/cloud_opacities/rf_data/ti_rad_1.txt', unpack=True)
	

extinct_ti = zeros((len(wavelength),nbin))
extinct_mg = zeros((len(wavelength),nbin))
extinct_fe = zeros((len(wavelength),nbin))
extinct_cr = zeros((len(wavelength),nbin))
extinct_mn = zeros((len(wavelength),nbin))
extinct_na = zeros((len(wavelength),nbin))
extinct_zn = zeros((len(wavelength),nbin))
extinct_k = zeros((len(wavelength),nbin))
extinct_al = zeros((len(wavelength),nbin))

iwave,ibin,ext_tio2, tio2_ssa, tio2_g =    np.loadtxt('/auto/home/dkpowell/EXO_TRANSMIT/UPDATED_Exotransmit/cloud_opacities/miex_data_morebins/tio2_opticalproperties.txt', unpack = True)
iwave,ibin,ext_mns, mns_ssa, mns_g =       np.loadtxt('/auto/home/dkpowell/EXO_TRANSMIT/UPDATED_Exotransmit/cloud_opacities/miex_data_morebins/mns_opticalproperties.txt', unpack = True)
iwave,ibin,ext_cr, cr_ssa, cr_g =          np.loadtxt('/auto/home/dkpowell/EXO_TRANSMIT/UPDATED_Exotransmit/cloud_opacities/miex_data_morebins/cr_opticalproperties.txt', unpack = True)
iwave,ibin,ext_al2o3, al2o3_ssa, al2o3_g = np.loadtxt('/auto/home/dkpowell/EXO_TRANSMIT/UPDATED_Exotransmit/cloud_opacities/miex_data_morebins/al2o3_opticalproperties.txt', unpack = True)
iwave,ibin,ext_na2s, na2s_ssa, na2s_g =    np.loadtxt('/auto/home/dkpowell/EXO_TRANSMIT/UPDATED_Exotransmit/cloud_opacities/miex_data_morebins/na2s_opticalproperties.txt', unpack = True)
iwave,ibin,ext_mg2sio4, mg2sio4_ssa, mg2sio4_g = np.loadtxt('/auto/home/dkpowell/EXO_TRANSMIT/UPDATED_Exotransmit/cloud_opacities/miex_data_morebins/mg2sio4_opticalproperties.txt', unpack = True)
iwave,ibin,ext_kcl, kcl_ssa, kcl_g =       np.loadtxt('/auto/home/dkpowell/EXO_TRANSMIT/UPDATED_Exotransmit/cloud_opacities/miex_data_morebins/kcl_opticalproperties.txt', unpack = True)
iwave,ibin,ext_zns, zns_ssa, zns_g =       np.loadtxt('/auto/home/dkpowell/EXO_TRANSMIT/UPDATED_Exotransmit/cloud_opacities/miex_data_morebins/zns_opticalproperties.txt', unpack = True)
iwave,ibin,ext_fe, fe_ssa, fe_g =          np.loadtxt('/auto/home/dkpowell/EXO_TRANSMIT/UPDATED_Exotransmit/cloud_opacities/miex_data_morebins/fe_opticalproperties.txt', unpack = True)

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
		
##########################################################################################
##########################################################################################
print 'Set it up!'
##########################################################################################
##########################################################################################
numx=1000
numz0=100
Rplanet = 6.99e7 #m
model_height = max(z)
z0_array = np.linspace(z[len(z)-2],z[0],num=numz0)

f_tio2 = interp2d(r[:,0],z,t_tio2)
f_het_fe = interp2d(r[:,2],z,t_het_fe) #nz,nbin
f_het_mg2sio4 = interp2d(r[:,3],z,t_het_mg2sio4)
f_pure_fe = interp2d(r[:,7],z,t_pure_fe)
f_cr = interp2d(r[:,8],z,t_cr)
f_het_cr = interp2d(r[:,8],z,t_het_cr)
f_mns = interp2d(r[:,5],z,t_mns)
f_na2s = interp2d(r[:,6],z,t_na2s)
f_zns = interp2d(r[:,10],z,t_zns)
f_kcl = interp2d(r[:,9],z,t_kcl)
f_al2o3 = interp2d(r[:,1],z,t_al2o3)

f_ext_ti = interp2d(r[:,0],wavelength,extinct_ti) #idw,nbin
f_ext_fe = interp2d(r[:,2],wavelength,extinct_fe)
f_ext_mg = interp2d(r[:,3],wavelength,extinct_mg)
f_ext_cr = interp2d(r[:,8],wavelength,extinct_cr)
f_ext_mn = interp2d(r[:,5],wavelength,extinct_mn)
f_ext_na = interp2d(r[:,6],wavelength,extinct_na)
f_ext_zn = interp2d(r[:,10],wavelength,extinct_zn)
f_ext_kc = interp2d(r[:,9],wavelength,extinct_k)
f_ext_al = interp2d(r[:,1],wavelength,extinct_al)

wavelength_sample = np.logspace(np.log10(min(wavelength)), np.log10(max(wavelength)),10)
	
		
od_ti1 = zeros((numz0,nbin))
od_mg1 = zeros((numz0,nbin))
od_fe1 = zeros((numz0,nbin))
od_cr1 = zeros((numz0,nbin))
od_mn1 = zeros((numz0,nbin))
od_na1 = zeros((numz0,nbin))
od_zn1 = zeros((numz0,nbin))
od_k1 = zeros((numz0,nbin))
od_al1 = zeros((numz0,nbin))

od_ti = zeros((numz0,len(wavelength_sample)))
od_mg = zeros((numz0,len(wavelength_sample)))
od_fe = zeros((numz0,len(wavelength_sample)))
od_cr = zeros((numz0,len(wavelength_sample)))
od_mn = zeros((numz0,len(wavelength_sample)))
od_na = zeros((numz0,len(wavelength_sample)))
od_zn = zeros((numz0,len(wavelength_sample)))
od_k = zeros((numz0,len(wavelength_sample)))
od_al = zeros((numz0,len(wavelength_sample)))

opt_depth = zeros((numz0,len(wavelength_sample)))

master_z = []
dx_array = []

for idx1,z0 in enumerate(z0_array):
	max_x = np.sqrt((Rplanet+model_height)**2.-(Rplanet+z0)**2.)
	dx = max_x/numx
	x_array = np.linspace(0,max_x,num=numx)
	z_array = np.sqrt((Rplanet+z0)**2.+x_array**2)-Rplanet
	master_z.append(z_array)
	dx_array.append(dx)
	
##########################################################################################
##########################################################################################	
print 'BIG LOOP TIME!'
##########################################################################################
##########################################################################################
op_tot = []
op_ti = []
op_mg = []
op_fe = []
op_cr = []
op_mn = []
op_na = []
op_zn = []
op_k = []
op_al = []
	
for idw,w in enumerate(wavelength_sample):
	for idx1,z_array in enumerate(master_z):
		for j in range(nbin):	
			od_ti1[idx1,j] = sum(f_tio2(r[j,0],z_array)*f_ext_ti(r[j,0],w)*100.*dx_array[idx1])
			od_mg1[idx1,j] = sum(f_het_mg2sio4(r[j,3],z_array)*f_ext_mg(r[j,3],w)*100.*dx_array[idx1])
			od_mn1[idx1,j] = sum(f_mns(r[j,5],z_array)*f_ext_mn(r[j,5],w)*100.*dx_array[idx1])
			od_na1[idx1,j] = sum(f_na2s(r[j,8],z_array)*f_ext_na(r[j,8],w)*100.*dx_array[idx1])
			od_zn1[idx1,j] = sum(f_zns(r[j,10],z_array)*f_ext_zn(r[j,10],w)*100.*dx_array[idx1])
			od_k1[idx1,j] = sum(f_kcl(r[j,9],z_array)*f_ext_kc(r[j,9],w)*100.*dx_array[idx1])
			od_al1[idx1,j] = sum(f_al2o3(r[j,1],z_array)*f_ext_al(r[j,1],w)*100.*dx_array[idx1])
			od_fe1[idx1,j] = sum(f_het_fe(r[j,2],z_array)*f_ext_fe(r[j,2],w)*100.*dx_array[idx1])
			od_cr1[idx1,j] = sum(f_het_cr(r[j,8],z_array)*f_ext_cr(r[j,8],w)*100.*dx_array[idx1])+sum(f_cr(r[j,8],z_array)*f_ext_cr(r[j,8],w)*100.*dx_array[idx1])
				
		
			
		od_ti[idx1,idw] = sum(od_ti1[idx1,:])
		od_mg[idx1,idw] = sum(od_mg1[idx1,:])
		od_mn[idx1,idw] = sum(od_mn1[idx1,:])
		od_na[idx1,idw] = sum(od_na1[idx1,:])
		od_zn[idx1,idw] = sum(od_zn1[idx1,:])
		od_k[idx1,idw] = sum(od_k1[idx1,:])
		od_al[idx1,idw] = sum(od_al1[idx1,:])
		od_fe[idx1,idw] = sum(od_fe1[idx1,:])
		od_cr[idx1,idw] = sum(od_cr1[idx1,:])
		opt_depth[idx1,idw] = sum(od_ti1[idx1,:]+od_mg1[idx1,:]+od_mn1[idx1,:]+od_na1[idx1,:]+od_zn1[idx1,:]+od_k1[idx1,:]
								+od_al1[idx1,:]+od_fe1[idx1,:]+od_cr1[idx1,:])
		
		
	if (od_fe[:,idw]>1.).any():
		blank = []
		zblank = []
		for i,o in enumerate(od_fe[:,idw]):
			blank.append(o)
			zblank.append(z0_array[i])
			if o > 1.:
				break
		f_opaque = interp1d(blank,zblank)
		op_fe.append(f_opaque(1.))
	else:
		op_fe.append(z[0])
			
	if (od_cr[:,idw]>1.).any():
		blank = []
		zblank = []
		for i,o in enumerate(od_cr[:,idw]):
			blank.append(o)
			zblank.append(z0_array[i])
			if o > 1.:
				break
		f_opaque = interp1d(blank,zblank)
		op_cr.append(f_opaque(1.))
	else:
		op_cr.append(z[0])
			
			
		
	if (od_ti[:,idw]>1.).any():
		blank = []
		zblank = []
		for i,o in enumerate(od_ti[:,idw]):
			blank.append(o)
			zblank.append(z0_array[i])
			if o > 1.:
				break
		f_opaque = interp1d(blank,zblank)
		op_ti.append(f_opaque(1.))
	else:
		op_ti.append(z[0])
			
	if (od_mg[:,idw]>1.).any():
		blank = []
		zblank = []
		for i,o in enumerate(od_mg[:,idw]):
			blank.append(o)
			zblank.append(z0_array[i])
			if o > 1.:
				break
		f_opaque = interp1d(blank,zblank)
		op_mg.append(f_opaque(1.))
	else:
		op_mg.append(z[0])
			
			
	if (od_mn[:,idw]>1.).any():
		blank = []
		zblank = []
		for i,o in enumerate(od_mn[:,idw]):
			blank.append(o)
			zblank.append(z0_array[i])
			if o > 1.:
				break
		f_opaque = interp1d(blank,zblank)
		op_mn.append(f_opaque(1.))
	else:
		op_mn.append(z[0])
			
	if (od_na[:,idw]>1.).any():
		blank = []
		zblank = []
		for i,o in enumerate(od_na[:,idw]):
			blank.append(o)
			zblank.append(z0_array[i])
			if o > 1.:
				break
		f_opaque = interp1d(blank,zblank)
		op_na.append(f_opaque(1.))
	else:
		op_na.append(z[0])
			
	if (od_zn[:,idw]>1.).any():
		blank = []
		zblank = []
		for i,o in enumerate(od_zn[:,idw]):
			blank.append(o)
			zblank.append(z0_array[i])
			if o > 1.:
				break
		f_opaque = interp1d(blank,zblank)
		op_zn.append(f_opaque(1.))
	else:
		op_zn.append(z[0])
			
	if (od_k[:,idw]>1.).any():
		blank = []
		zblank = []
		for i,o in enumerate(od_k[:,idw]):
			blank.append(o)
			zblank.append(z0_array[i])
			if o > 1.:
				break
		f_opaque = interp1d(blank,zblank)
		op_k.append(f_opaque(1.))
	else:
		op_k.append(z[0])
			
	if (od_al[:,idw]>1.).any():
		blank = []
		zblank = []
		for i,o in enumerate(od_al[:,idw]):
			blank.append(o)
			zblank.append(z0_array[i])
			if o > 1.:
				break
		f_opaque = interp1d(blank,zblank)
		op_al.append(f_opaque(1.))
	else:
		op_al.append(z[0])
		
		
		
	if (opt_depth[:,idw]>1.).any():
		blank = []
		zblank = []
		for i,o in enumerate(opt_depth[:,idw]):
			blank.append(o)
			zblank.append(z0_array[i])
			if o > 1.:
				break
		f_opaque = interp1d(blank,zblank)
		op_tot.append(f_opaque(1.))
	else:
		op_tot.append(z[0])


f_pressure = interp1d(z,np.log10(p))
tot = []
ti = []
mg = []
fe = []
cr = []
mn = []
na = []
kcl = []
al = []
zn = []


for zed in op_tot:
	tot.append(10.**f_pressure(zed))
tot = np.array(tot)/1.e6
for zed in op_ti:
	ti.append(10.**f_pressure(zed))
ti = np.array(ti)/1.e6
for zed in op_mg:
	mg.append(10.**f_pressure(zed))
mg = np.array(mg)/1.e6
for zed in op_fe:
	fe.append(10.**f_pressure(zed))
fe = np.array(fe)/1.e6
for zed in op_cr:
	cr.append(10.**f_pressure(zed))
cr = np.array(cr)/1.e6
for zed in op_mn:
	mn.append(10.**f_pressure(zed))
mn = np.array(mn)/1.e6
for zed in op_na:
	na.append(10.**f_pressure(zed))
na = np.array(na)/1.e6
for zed in op_k:
	kcl.append(10.**f_pressure(zed))
kcl = np.array(kcl)/1.e6
for zed in op_al:
	al.append(10.**f_pressure(zed))
al = np.array(al)/1.e6
for zed in op_zn:
	zn.append(10.**f_pressure(zed))
zn = np.array(zn)/1.e6
	

print tot




np.savetxt(path+'/'+name+'_opticaldepth_level.txt',np.c_[wavelength_sample,tot])


		
			