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
from scipy.interpolate import interp1d, interp2d, RectBivariateSpline

# Read in inputs from command line
opts, args = getopt.getopt(sys.argv[1:], "h")
file = args[0]
name = args[1]
path = args[2]

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
time = zeros(nstep/iskip)

tio2, het_fe, het_mg2sio4, pure_fe, pure_mg2sio4, cr, mns, na2s,zns,kcl,al2o3 = zeros((nz,nbin,nstep/iskip)),zeros((nz,nbin,nstep/iskip)), zeros((nz,nbin,nstep/iskip)), zeros((nz,nbin,nstep/iskip)), zeros((nz,nbin,nstep/iskip)), zeros((nz,nbin,nstep/iskip)), zeros((nz,nbin,nstep/iskip)), zeros((nz,nbin,nstep/iskip)), zeros((nz,nbin,nstep/iskip)), zeros((nz,nbin,nstep/iskip)), zeros((nz,nbin,nstep/iskip))
t_tio2, t_het_fe, t_het_mg2sio4, t_pure_fe, t_pure_mg2sio4, t_cr, t_mns, t_na2s,t_zns,t_kcl,t_al2o3 = zeros((nz,nbin)),zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin))
count = 0

M_KCL = 1.98
M_ZNS = 4.1
M_NA2S = 1.856
M_MNS = 3.3
M_CR = 7.19
M_FE = 7.874
M_MG2SIO4 = 3.21
M_TIO2 = 4.23
M_AL2O3 = 3.95

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
	
			
for i in range(nstep/iskip-1):
	time[i] = float(data.next())
	for j in range(nbin):
		for k in range(nz):
			line = data.next().split()
			
			tio2[k,j,i] = float(line[2])
			
			het_fe[k,j,i] = float(line[3])
			
			het_mg2sio4[k,j,i] = float(line[4])
			
			pure_fe[k,j,i] = float(line[5])
			
			pure_mg2sio4[k,j,i] = float(line[6])
			
			cr[k,j,i] = float(line[7])
			
			mns[k,j,i] = float(line[8])
			
			na2s[k,j,i] = float(line[9])
			
			zns[k,j,i] = float(line[10])
			
			kcl[k,j,i] = float(line[11])
			
			al2o3[k,j,i] = float(line[12])
			
	if i > nstep/iskip-10 :
		count += 1
		
		t_tio2[:,:] += tio2[:,:,i]
		
		t_het_fe[:,:] += het_fe[:,:,i]
		
		t_het_mg2sio4[:,:] += het_mg2sio4[:,:,i]
		
		t_pure_fe[:,:] += pure_fe[:,:,i]
		
		t_pure_mg2sio4[:,:] += pure_mg2sio4[:,:,i]
		
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
t_pure_mg2sio4 /= count
t_cr /= count
t_mns /= count
t_na2s /= count
t_zns /= count
t_kcl /= count
t_al2o3 /= count


wavelength,extinct_ti = np.loadtxt('../../EXO_TRANSMIT/Exo_Transmit-master/cloud_opacities/rf_data/ti_rad_1.txt', unpack=True)		

od_ti = zeros((nz,nbin))
od_mg = zeros((nz,nbin))
od_fe = zeros((nz,nbin))
od_cr = zeros((nz,nbin))
od_mn = zeros((nz,nbin))
od_na = zeros((nz,nbin))
od_zn = zeros((nz,nbin))
od_k = zeros((nz,nbin))
od_al = zeros((nz,nbin))

opt_depth = zeros((nz, len(wavelength)))

extinct_ti = zeros((len(wavelength),nbin))
extinct_mg = zeros((len(wavelength),nbin))
extinct_fe = zeros((len(wavelength),nbin))
extinct_cr = zeros((len(wavelength),nbin))
extinct_mn = zeros((len(wavelength),nbin))
extinct_na = zeros((len(wavelength),nbin))
extinct_zn = zeros((len(wavelength),nbin))
extinct_k = zeros((len(wavelength),nbin))
extinct_al = zeros((len(wavelength),nbin))
	
for j in range(nbin):
	wi,exti = np.loadtxt('../../EXO_TRANSMIT/Exo_Transmit-master/cloud_opacities/rf_data/ti_rad_'+str(j)+'.txt', unpack=True)
	wi,exmg = np.loadtxt('../../EXO_TRANSMIT/Exo_Transmit-master/cloud_opacities/rf_data/mg_rad_'+str(j)+'.txt', unpack=True)
	wi,exfe = np.loadtxt('../../EXO_TRANSMIT/Exo_Transmit-master/cloud_opacities/rf_data/fe_rad_'+str(j)+'.txt', unpack=True)
	wi,excr = np.loadtxt('../../EXO_TRANSMIT/Exo_Transmit-master/cloud_opacities/rf_data/cr_rad_'+str(j)+'.txt', unpack=True)
	wi,exmn = np.loadtxt('../../EXO_TRANSMIT/Exo_Transmit-master/cloud_opacities/rf_data/mn_rad_'+str(j)+'.txt', unpack=True)
	wi,exna = np.loadtxt('../../EXO_TRANSMIT/Exo_Transmit-master/cloud_opacities/rf_data/na_rad_'+str(j)+'.txt', unpack=True)
	wi,exzn = np.loadtxt('../../EXO_TRANSMIT/Exo_Transmit-master/cloud_opacities/rf_data/zn_rad_'+str(j)+'.txt', unpack=True)
	wi,exk = np.loadtxt('../../EXO_TRANSMIT/Exo_Transmit-master/cloud_opacities/rf_data/k_rad_'+str(j)+'.txt', unpack=True)
	wi,exal = np.loadtxt('../../EXO_TRANSMIT/Exo_Transmit-master/cloud_opacities/rf_data/al_rad_'+str(j)+'.txt', unpack=True)
	extinct_ti[:,j] = exti
	extinct_mg[:,j] = exmg
	extinct_fe[:,j] = exfe
	extinct_cr[:,j] = excr
	extinct_mn[:,j] = exmn
	extinct_na[:,j] = exna
	extinct_k[:,j] = exk
	extinct_al[:,j] = exal
	extinct_zn[:,j] = exzn

opaque = []		
		
for idw,w in enumerate(wavelength):	
	for k in range(nz):
		for j in range(nbin):	
			dz = (z[k]-z[k-1])*100.
			od_ti[k,j] = t_tio2[k,j]*extinct_ti[idw,j]*dz
			od_mg[k,j] = (t_het_mg2sio4[k,j]*extinct_mg[idw,j]*dz)+(t_pure_mg2sio4[k,j]*extinct_mg[idw,j]*dz)
			od_fe[k,j] = (t_het_fe[k,j]*extinct_fe[idw,j]*dz)+(t_pure_fe[k,j]*extinct_fe[idw,j]*dz)
			od_cr[k,j] = t_cr[k,j]*extinct_cr[idw,j]*dz
			od_mn[k,j] = t_mns[k,j]*extinct_mn[idw,j]*dz
			od_na[k,j] = t_na2s[k,j]*extinct_na[idw,j]*dz
			od_k[k,j] = t_kcl[k,j]*extinct_k[idw,j]*dz
			od_al[k,j] = t_al2o3[k,j]*extinct_al[idw,j]*dz
			od_zn[k,j] = t_zns[k,j]*extinct_zn[idw,j]*dz
		
		opt_depth[k,idw] = sum(od_ti[k,:]+od_mg[k,:]+od_fe[k,:]+od_cr[k,:] + od_mn[k,:] + od_na[k,:] + od_k[k,:] + od_al[k,:] + od_zn[k,:])
	opt_depth[:,idw] = np.cumsum(opt_depth[:,idw][::-1])
	opt_depth[:,idw] = opt_depth[:,idw][::-1]
	
	if (opt_depth[:,idw]>1.).any():
		blank = []
		pblank = []
		p_inver = p[::-1]
		for i,o in enumerate(opt_depth[:,idw][::-1]):
			blank.append(o)
			pblank.append(p_inver[i])
			if o > 1.:
				break
		f_opaque = interp1d(blank,pblank)
		opaque.append(f_opaque(1.)/1.e6)
	else:
		opaque.append(z[0])		

combined = np.vstack((wavelength,opaque)).T	
np.savetxt(path+'/ocl_'+name+'.txt', combined)	

f_ocl = interp1d(wavelength,opaque)
print f_ocl(1.25)

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
		
X,Y = np.meshgrid(wavelength, p/1.e6)
f,ax = plt.subplots(1,1)
ax.set_axis_bgcolor('k')
im = ax.pcolormesh(X,Y,opt_depth, cmap=cmaps.magma, norm=LogNorm())
im.set_clim(vmin=1.e-7,vmax=1.e1)
cbar = plt.colorbar(im)
ax.scatter(wavelength,opaque, c='w')
ax.plot(wavelength,opaque, c='w')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim(max(p)/1.e6,min(p)/1.e6)
ax.set_xlim(min(wavelength),max(wavelength))
ax.set_ylabel(r'Pressure [bar]')
ax.set_xlabel(r'Wavelength [$\mu$m]')
cbar.set_label(r'$\tau$', rotation=-90, labelpad=15)
plt.savefig(path+'/nadir_spec_'+name+'.jpg', format = 'jpg', bbox_inches='tight')
plt.close()
		
		
		
		
		
		
		
		
		