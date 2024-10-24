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
mass_tio2, mass_het_fe, mass_het_mg2sio4, mass_pure_fe, mass_pure_mg2sio4, mass_cr, mass_mns, mass_na2s, mass_zns, mass_kcl, mass_al2o3 = zeros((nz,nbin,nstep/iskip)),zeros((nz,nbin,nstep/iskip)), zeros((nz,nbin,nstep/iskip)), zeros((nz,nbin,nstep/iskip)), zeros((nz,nbin,nstep/iskip)), zeros((nz,nbin,nstep/iskip)), zeros((nz,nbin,nstep/iskip)), zeros((nz,nbin,nstep/iskip)), zeros((nz,nbin,nstep/iskip)), zeros((nz,nbin,nstep/iskip)), zeros((nz,nbin,nstep/iskip))

t_tio2, t_het_fe, t_het_mg2sio4, t_pure_fe, t_pure_mg2sio4, t_cr, t_mns, t_na2s,t_zns,t_kcl,t_al2o3 = zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin))
t_mass_tio2, t_mass_het_fe, t_mass_het_mg2sio4, t_mass_pure_fe, t_mass_pure_mg2sio4, t_mass_cr, t_mass_mns, t_mass_na2s, t_mass_zns, t_mass_kcl, t_mass_al2o3 = zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin))

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
			
for i in range(nstep/iskip-1):
	time[i] = float(data.next())
	for j in range(nbin):
		for k in range(nz):
			line = data.next().split()
			
			tio2[k,j,i] = float(line[2])/(math.log(ru[j,1]) - math.log(rl[j,1]))
			mass_tio2[k,j,i] = tio2[k,j,i]*(M_TIO2*((4./3.)*np.pi*(r[j,1]*1e-4)**3.))
			
			het_fe[k,j,i] = float(line[3])/(math.log(ru[j,1]) - math.log(rl[j,1]))
			mass_het_fe[k,j,i] = het_fe[k,j,i]*(M_FE*((4./3.)*np.pi*(r[j,1]*1e-4)**3.))
			
			het_mg2sio4[k,j,i] = float(line[4])/(math.log(ru[j,1]) - math.log(rl[j,1]))
			mass_het_mg2sio4[k,j,i] = het_mg2sio4[k,j,i]*(M_MG2SIO4*((4./3.)*np.pi*(r[j,1]*1e-4)**3.))
			
			pure_fe[k,j,i] = float(line[5])/(math.log(ru[j,1]) - math.log(rl[j,1]))
			mass_pure_fe[k,j,i] = pure_fe[k,j,i]*(M_FE*((4./3.)*np.pi*(r[j,1]*1e-4)**3.))
			
			pure_mg2sio4[k,j,i] = float(line[6])/(math.log(ru[j,1]) - math.log(rl[j,1]))
			mass_pure_mg2sio4[k,j,i] = pure_mg2sio4[k,j,i]*(M_MG2SIO4*((4./3.)*np.pi*(r[j,1]*1e-4)**3.))
			
			cr[k,j,i] = float(line[7])/(math.log(ru[j,1]) - math.log(rl[j,1]))
			mass_cr[k,j,i] = cr[k,j,i]*(M_CR*((4./3.)*np.pi*(r[j,1]*1e-4)**3.))
			
			mns[k,j,i] = float(line[8])/(math.log(ru[j,1]) - math.log(rl[j,1]))
			mass_mns[k,j,i] = mns[k,j,i]*(M_MNS*((4./3.)*np.pi*(r[j,1]*1e-4)**3.))
			
			na2s[k,j,i] = float(line[9])/(math.log(ru[j,1]) - math.log(rl[j,1]))
			mass_na2s[k,j,i] = na2s[k,j,i]*(M_NA2S*((4./3.)*np.pi*(r[j,1]*1e-4)**3.))
			
			zns[k,j,i] = float(line[10])/(math.log(ru[j,1]) - math.log(rl[j,1]))
			mass_zns[k,j,i] = zns[k,j,i]*(M_ZNS*((4./3.)*np.pi*(r[j,1]*1e-4)**3.))
			
			kcl[k,j,i] = float(line[11])/(math.log(ru[j,1]) - math.log(rl[j,1]))
			mass_kcl[k,j,i] = kcl[k,j,i]*(M_KCL*((4./3.)*np.pi*(r[j,1]*1e-4)**3.))
			
			al2o3[k,j,i] = float(line[12])/(math.log(ru[j,1]) - math.log(rl[j,1]))
			mass_al2o3[k,j,i] = al2o3[k,j,i]*(M_AL2O3*((4./3.)*np.pi*(r[j,1]*1e-4)**3.))
			
	if i > nstep/iskip - 10:
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
		
		t_mass_tio2[:,:] += mass_tio2[:,:,i]
		t_mass_het_fe[:,:] += mass_het_fe[:,:,i]
		t_mass_het_mg2sio4[:,:] += mass_het_mg2sio4[:,:,i]
		t_mass_pure_fe[:,:] += mass_pure_fe[:,:,i]
		t_mass_pure_mg2sio4[:,:] += mass_pure_mg2sio4[:,:,i]
		t_mass_cr[:,:] += mass_cr[:,:,i]
		t_mass_mns[:,:] += mass_mns[:,:,i]
		t_mass_na2s[:,:] += mass_na2s[:,:,i]
		t_mass_zns[:,:] += mass_zns[:,:,i]
		t_mass_kcl[:,:] += mass_kcl[:,:,i]
		t_mass_al2o3[:,:] += mass_al2o3[:,:,i]
				

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

t_mass_tio2 /= count
t_mass_het_fe /= count
t_mass_het_mg2sio4 /= count
t_mass_pure_fe /= count
t_mass_pure_mg2sio4 /= count
t_mass_cr /= count
t_mass_mns /= count
t_mass_na2s /= count
t_mass_zns /= count
t_mass_kcl /= count
t_mass_al2o3 /= count

tot_num = t_tio2+t_het_fe+t_het_mg2sio4+t_pure_fe+t_pure_mg2sio4+t_cr+t_mns+t_na2s+t_zns+t_kcl+t_al2o3

tot_mass = t_mass_tio2+t_mass_het_fe+t_mass_het_mg2sio4+t_mass_pure_fe+t_mass_pure_mg2sio4+t_mass_cr+t_mass_mns+t_mass_na2s+t_mass_zns+t_mass_kcl+t_mass_al2o3

##########################################################################################
##########################################################################################
print 'Plotting it up good!'
##########################################################################################
##########################################################################################

print 'Making Figure ', 1
plt.rcParams['font.size'] = 16

f,ax = plt.subplots(1,1)
ax.set_axis_bgcolor('k')	

X, Y = np.meshgrid(r[:,0],p[:]/1.e6)

im = ax.pcolormesh(X,Y,tot_num, norm=LogNorm(), cmap=cmaps.magma)
im.set_clim(vmin=1.e-5,vmax=1.e5)

cbar = plt.colorbar(im)
cbar.set_label(r'dN/dLn(r) [cm$^{-3}$]', fontsize=16, rotation=-90, labelpad=30)

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim(1.e2,min(p)/1.e6)
ax.set_xlim(min(r[:,0]),max(r[:,0]))
ax.set_ylabel(r'Pressure [bar]')
ax.set_xlabel(r'Particle Radius [$\mu$m]')
plt.axhline(0.00341895869576, c='w', lw=2, ls='--')

savefig(path+'/'+name+'_tot_numvert.jpg', format = 'jpg', bbox_inches='tight')
plt.close()
	
print 'Making Figure ', 2
plt.rcParams['font.size'] = 16

f,ax = plt.subplots(1,1)
ax.set_axis_bgcolor('k')	

X, Y = np.meshgrid(r[:,0],p[:]/1.e6)

im = ax.pcolormesh(X,Y,tot_mass, norm=LogNorm(), cmap=cmaps.magma)
im.set_clim(vmin=1.e-18,vmax=1.e-9)

cbar = plt.colorbar(im)
cbar.set_label(r'dM/dLn(r) [g cm$^{-3}$]', fontsize=16, rotation=-90, labelpad=30)

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim(1.e2,min(p)/1.e6)
ax.set_xlim(min(r[:,0]),max(r[:,0]))
ax.set_ylabel(r'Pressure [bar]')
ax.set_xlabel(r'Particle Radius [$\mu$m]')
plt.axhline(0.00341895869576, c='w', lw=2, ls='--')

savefig(path+'/'+name+'_tot_massvert.jpg', format = 'jpg', bbox_inches='tight')
plt.close()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	