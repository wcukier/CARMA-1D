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
stprnt = float(args[3])
pskip = int(args[4])

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

mix_ti, svp_ti, mix_fe, svp_fe, mix_mg, svp_mg, mix_cr, svp_cr, mix_al, svp_al = zeros((nz,ntime)), zeros((nz,ntime)), zeros((nz,ntime)), zeros((nz,ntime)), zeros((nz,ntime)), zeros((nz,ntime)), zeros((nz,ntime)), zeros((nz,ntime)), zeros((nz,ntime)), zeros((nz,ntime))

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
			
			mix_ti[k,i] = float(line[13])
			svp_ti[k,i] = float(line[14])
			mix_fe[k,i] = float(line[15])
			svp_fe[k,i] = float(line[16])
			mix_mg[k,i] = float(line[17])
			svp_mg[k,i] = float(line[18])
			mix_cr[k,i] = float(line[19])
			svp_cr[k,i] = float(line[20])
			mix_al[k,i] = float(line[29])
			svp_al[k,i] = float(line[30])
	
	
	if i >= stprnt and i % pskip == 0:
	
		print('Making Figure '+str(i+1)+' out of '+str(nstep/iskip))
		plt.rcParams['figure.figsize'] = (6,5)
		plt.rcParams['legend.frameon'] = False
		plt.rcParams['legend.fontsize'] = 14
		plt.rcParams['legend.borderpad'] = 0.1
		plt.rcParams['legend.labelspacing'] = 0.1
		plt.rcParams['legend.handletextpad'] = 0.1
		plt.rcParams['font.family'] = 'stixgeneral'
		plt.rcParams['font.size'] = 16
		
	
		fig, ax = subplots(nrows = 1, ncols = 1)
		ax.plot(mix_ti[:,i],z, c='b')
		ax.plot(svp_ti[:,i],z, c='b', ls='--')
		
		ax.plot(mix_fe[:,i],z, c='r')
		ax.plot(svp_fe[:,i],z, c='r', ls='--')
		
		ax.set_yscale('log')
		ax.set_xscale('log')
		#ax.set_ylim(1.e2,min(p)/1.e6)
		#ax.set_xlim(min(r[:,0]),max(r[:,0]))
		ax.set_ylabel(r'Height')
		ax.set_xlabel(r'gas mmr')
		plt.savefig(path+'/gas_mmr_'+str(name)+'_'+str(i)+'.png', format = 'png', bbox_inches='tight')
		plt.close()
	
		print('Making Figure '+str(i+1)+' out of '+str(nstep/iskip))
		plt.rcParams['figure.figsize'] = (6,5)
		plt.rcParams['legend.frameon'] = False
		plt.rcParams['legend.fontsize'] = 14
		plt.rcParams['legend.borderpad'] = 0.1
		plt.rcParams['legend.labelspacing'] = 0.1
		plt.rcParams['legend.handletextpad'] = 0.1
		plt.rcParams['font.family'] = 'stixgeneral'
		plt.rcParams['font.size'] = 16
		
	
		fig, ax = subplots(nrows = 1, ncols = 1)
		ax.plot(mix_ti[:,i],p/1.e6, c='b')
		ax.plot(svp_ti[:,i],p/1.e6, c='b', ls='--')
		
		ax.plot(mix_fe[:,i],p/1.e6, c='r')
		ax.plot(svp_fe[:,i],p/1.e6, c='r', ls='--')
		
		ax.set_yscale('log')
		ax.set_xscale('log')
		ax.set_ylim(1.e2,min(p)/1.e6)
		#ax.set_xlim(min(r[:,0]),max(r[:,0]))
		ax.set_ylabel(r'Height')
		ax.set_xlabel(r'gas mmr')
		plt.savefig(path+'/gas_mmr_pressure_'+str(name)+'_'+str(i)+'.png', format = 'png', bbox_inches='tight')
		plt.close()


infile.close()









