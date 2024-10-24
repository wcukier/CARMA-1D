#! /usr/bin/env python3.6
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

#./cloud_layers.py ../run/carma/hjt_west_1800.txt 1800_west layer_plots/ > log.txt &

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
print('Read in output file')
infile=open(file,'r')
line = infile.readline().split()

nz,ngroup,nelem,nbin,ngas,nstep,iskip = map(int,line)
r,ms,dr,rl,ru = zeros((nbin,ngroup)),zeros((nbin,ngroup)),zeros((nbin,ngroup)),zeros((nbin,ngroup)),zeros((nbin,ngroup))
p,t,z = zeros(nz),zeros(nz),zeros(nz)

ntime=int(nstep/iskip)

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
			
#np.savetxt('new_ti_rad.txt', r[:,0].T)
#np.savetxt('new_al_rad.txt', r[:,1].T)
#np.savetxt('new_fe_rad.txt', r[:,2].T)
#np.savetxt('new_mg_rad.txt', r[:,3].T)
#np.savetxt('new_cr_rad.txt', r[:,4].T)
#np.savetxt('new_mn_rad.txt', r[:,5].T)
#np.savetxt('new_na_rad.txt', r[:,6].T)
#np.savetxt('new_pfe_rad.txt', r[:,7].T)
#np.savetxt('new_pcr_rad.txt', r[:,8].T)
#np.savetxt('new_k_rad.txt', r[:,9].T)
#np.savetxt('new_zn_rad.txt', r[:,10].T)


for i in range(ntime):
	time[i] = float(infile.readline())
	#if i== 0:
	#	time[i] = float(infile.readline())
	for j in range(nbin):
		for k in range(nz):
			
			line = infile.readline().split()
			
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
			
#print max(map(max,pm_al)), max(map(max,pm_mg)), max(map(max,pm_fe+pm_pure_fe)), max(map(max,pm_cr+pm_pure_cr)), max(map(max,pm_ti)), max(map(max,pm_mn))
#print max(map(max,p_al)), max(map(max,p_mg)), max(map(max,p_fe+p_pure_fe)), max(map(max,p_cr+p_pure_cr)), max(map(max,p_ti)), max(map(max,p_mn))

#print max(map(max,t_al2o3)), max(map(max,t_het_mg2sio4)), max(map(max,t_het_fe+t_pure_fe)), max(map(max,t_cr+t_het_cr)), max(map(max,t_tio2)), max(map(max,t_mns))
#print max(map(max,t_mass_al2o3)), max(map(max,t_mass_het_mg2sio4)), max(map(max,t_mass_het_fe+t_mass_pure_fe)), max(map(max,t_mass_cr+t_mass_het_cr)), max(map(max,t_mass_tio2)), max(map(max,t_mass_mns))

X1, Y1 = np.meshgrid(r[:,0],p[:]/1.e6)
X2, Y2 = np.meshgrid(r[:,1],p[:]/1.e6)
X3, Y3 = np.meshgrid(r[:,2],p[:]/1.e6)
X4, Y4 = np.meshgrid(r[:,3],p[:]/1.e6)
X5, Y5 = np.meshgrid(r[:,4],p[:]/1.e6)
X6, Y6 = np.meshgrid(r[:,5],p[:]/1.e6)
X7, Y7 = np.meshgrid(r[:,6],p[:]/1.e6)
X8, Y8 = np.meshgrid(r[:,7],p[:]/1.e6)
X9, Y9 = np.meshgrid(r[:,8],p[:]/1.e6)
X10, Y10 = np.meshgrid(r[:,9],p[:]/1.e6)
X11, Y11 = np.meshgrid(r[:,10],p[:]/1.e6)

plt.rcParams['figure.figsize'] = (6,5)
plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['legend.borderpad'] = 0.1
plt.rcParams['legend.labelspacing'] = 0.1
plt.rcParams['legend.handletextpad'] = 0.1
plt.rcParams['font.family'] = 'stixgeneral'
plt.rcParams['font.size'] = 16

levels_num = [1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4]

fig, ax = subplots(nrows = 1, ncols = 1)
im1 = ax.contourf(X1,Y1,p_ti, norm=LogNorm(),levels=levels_num, cmap='Blues', alpha=0.5)
im2 = ax.contourf(X2,Y2,p_al, norm=LogNorm(),levels=levels_num, cmap='Greens', alpha=0.5)
im4 = ax.contourf(X4,Y4,p_mg, norm=LogNorm(),levels=levels_num, cmap='Purples', alpha=0.5)
im3 = ax.contourf(X3,Y3,p_fe, norm=LogNorm(),levels=levels_num, cmap='Reds', alpha=0.75)
im5 = ax.contourf(X5,Y5,p_cr+p_pure_cr, norm=LogNorm(),levels=levels_num, cmap='Oranges', alpha=0.75)
cbar = plt.colorbar(im1)
cbar.set_label(r'dN/dLn(r) [cm$^{-3}$]', fontsize=16, rotation=-90, labelpad=30)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim(1.e2,min(p)/1.e6)
ax.set_xlim(min(r[:,0]),max(r[:,0]))
ax.set_ylabel(r'Pressure [bar]')
ax.set_xlabel(r'Particle Radius [$\mu$m]')
plt.savefig(path+'/layers_num_'+str(name)+'.jpg', format = 'jpg', bbox_inches='tight')
plt.savefig(path+'/layers_num_'+str(name)+'.eps', format = 'eps', bbox_inches='tight')
plt.close()

levels_mass = [1e-14,1e-13,1e-12,1e-11,1e-10,1e-9]

fig, ax = subplots(nrows = 1, ncols = 1)
im1 = ax.contourf(X1,Y1,pm_ti, norm=LogNorm(),levels=levels_mass, cmap='Blues', alpha=0.5)
im2 = ax.contourf(X2,Y2,pm_al, norm=LogNorm(),levels=levels_mass, cmap='Greens', alpha=0.5)
im3 = ax.contourf(X3,Y3,pm_fe, norm=LogNorm(),levels=levels_mass, cmap='Reds', alpha=0.5)
im5 = ax.contourf(X4,Y4,pm_mg, norm=LogNorm(),levels=levels_mass, cmap='Purples', alpha=0.5)
im6 = ax.contourf(X5,Y5,pm_cr+pm_pure_cr, norm=LogNorm(),levels=levels_mass, cmap='Oranges', alpha=0.5)
cbar = plt.colorbar(im1)
cbar.set_label(r'dM/dLn(r) [g cm$^{-3}$]', fontsize=16, rotation=-90, labelpad=30)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim(1.e2,min(p)/1.e6)
ax.set_xlim(min(r[:,0]),max(r[:,0]))
ax.set_ylabel(r'Pressure [bar]')
ax.set_xlabel(r'Particle Radius [$\mu$m]')
plt.savefig(path+'/layers_mass_'+str(name)+'.jpg', format = 'jpg', bbox_inches='tight')
plt.savefig(path+'/layers_mass_'+str(name)+'.eps', format = 'eps', bbox_inches='tight')
plt.close()


infile.close()









