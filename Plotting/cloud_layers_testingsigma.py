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
			


for i in range(ntime):
	time[i] = float(infile.readline())
	#accidentally forgot to print out the time after the first step in the cartesian run
	#if i== 0:
	#	time[i] = float(infile.readline())
	
	for j in range(nbin):
		for k in range(nz):
			
			line = infile.readline().split()
			
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

levels_num = [1e-12,1e-10,1e-8,1e-6,1e-4]

fig, ax = subplots(nrows = 1, ncols = 1)
im1 = ax.contourf(X1,Y1,t_tio2, norm=LogNorm(),levels=levels_num, cmap='Blues', alpha=0.5)
im2 = ax.contourf(X2,Y2,t_al2o3, norm=LogNorm(),levels=levels_num, cmap='Greens', alpha=0.5)
im4 = ax.contourf(X4,Y4,t_het_mg2sio4, norm=LogNorm(),levels=levels_num, cmap='Purples', alpha=0.5)
im3 = ax.contourf(X3,Y3,t_het_fe+t_pure_fe, norm=LogNorm(),levels=levels_num, cmap='Reds', alpha=0.75)
im5 = ax.contourf(X5,Y5,t_het_cr+t_cr, norm=LogNorm(),levels=levels_num, cmap='Oranges', alpha=0.75)
cbar = plt.colorbar(im1)
cbar.set_label(r'MMR', fontsize=16, rotation=-90, labelpad=30)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim(1.e2,min(p)/1.e6)
ax.set_xlim(min(r[:,0]),max(r[:,0]))
ax.set_ylabel(r'Pressure [bar]')
ax.set_xlabel(r'Particle Radius [$\mu$m]')
plt.savefig(path+'/layers_mmr_'+str(name)+'.png', format = 'png', bbox_inches='tight')
plt.close()



infile.close()









