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

#./cloud_opacity_out.py hjt_west_1800.txt ../T_P/t_p_1800_west_limb.dat 1800_west opacities/ > log.txt &

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

'''
if name == '2000_boteast':
	ntime = 299
if name == '2000_botwest':
	ntime = 322
if name == '2000_topeast':
	ntime = 319
if name == '2000_topwest':
	ntime = 301
if name == '1800_east':
	ntime = 249
if name == '1900_east':
	ntime = 303
if name == '2000_east':
	ntime = 238
if name == '2100_east':
	ntime = 302
if name == '1800_west':
	ntime = 248
if name == '1900_west':
	ntime = 312
if name == '2000_west':
	ntime = 258
if name == '2100_west':
	ntime = 314
if name == '1800_pole':
	ntime = 249
if name == '1900_pole':
	ntime = 317
if name == '2000_pole':
	ntime = 249
if name == '2100_pole':
	ntime = 318
'''

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
			
			
	if i > ntime-30 :
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
	

extinct_ti = zeros((len(wavelength),nbin))
extinct_mg = zeros((len(wavelength),nbin))
extinct_fe = zeros((len(wavelength),nbin))
extinct_cr = zeros((len(wavelength),nbin))
extinct_mn = zeros((len(wavelength),nbin))
extinct_na = zeros((len(wavelength),nbin))
extinct_zn = zeros((len(wavelength),nbin))
extinct_k = zeros((len(wavelength),nbin))
extinct_al = zeros((len(wavelength),nbin))

iwave,ibin,ext_tio2, tio2_ssa, tio2_g =    np.loadtxt('miex_data/tio2_opticalproperties.txt', unpack = True)
iwave,ibin,ext_mns, mns_ssa, mns_g =       np.loadtxt('miex_data/mns_opticalproperties.txt', unpack = True)
iwave,ibin,ext_cr, cr_ssa, cr_g =          np.loadtxt('miex_data/cr_opticalproperties.txt', unpack = True)
iwave,ibin,ext_al2o3, al2o3_ssa, al2o3_g = np.loadtxt('miex_data/al2o3_opticalproperties.txt', unpack = True)
iwave,ibin,ext_na2s, na2s_ssa, na2s_g =    np.loadtxt('miex_data/na2s_opticalproperties.txt', unpack = True)
iwave,ibin,ext_mg2sio4, mg2sio4_ssa, mg2sio4_g = np.loadtxt('miex_data/mg2sio4_opticalproperties.txt', unpack = True)
iwave,ibin,ext_kcl, kcl_ssa, kcl_g =       np.loadtxt('miex_data/kcl_opticalproperties.txt', unpack = True)
iwave,ibin,ext_zns, zns_ssa, zns_g =       np.loadtxt('miex_data/zns_opticalproperties.txt', unpack = True)
iwave,ibin,ext_fe, fe_ssa, fe_g =          np.loadtxt('miex_data/fe_opticalproperties.txt', unpack = True)

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

plotty = 1

if plotty:
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



		
			