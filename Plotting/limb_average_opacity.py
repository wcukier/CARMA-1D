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

#./cloud_opacity_out_morebins.py ../../../New_ExoCARMA/apr7_half_latentheat/run/carma/hjt_west_2000.txt ../T_P/t_p_2000_west_limb.dat 2000_west apr11_halflatent_opacities/ > log.txt &
#./cloud_opacity_out_morebins.py ../../../New_ExoCARMA/apr7_half_latentheat/run/carma/hjt_east_2000.txt ../T_P/t_p_2000_east_limb.dat 2000_east apr11_halflatent_opacities/ > log.txt &
#./cloud_opacity_out_morebins.py ../../../New_ExoCARMA/apr7_half_latentheat/run/carma/hjt_pole_2000.txt ../T_P/t_p_2000_pole.dat 2000_pole apr11_halflatent_opacities/ > log.txt &

# Read in inputs from command line
opts, args = getopt.getopt(sys.argv[1:], "h")
file = args[0]
file1 = args[0]
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
print 'BIG LOOP TIME!'
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



f_tot_opaque = interp2d(wavelength,p/1.e6,tot_cloud_opaque)
i,pressure,temperature = np.loadtxt(tp_file, unpack=True, skiprows = 1)
pressure*=1e-5

print p/1.e6, pressure

output = f_tot_opaque(wavelength,pressure)
output *= 100. #opacity is in units of m^-1

np.savetxt(path+'/cloud_opacities_'+name+'.txt', output.T, fmt='%.3e')



		
			