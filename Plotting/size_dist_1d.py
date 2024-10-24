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
from scipy.stats import lognorm

# ./size_dist_1d.py ../run/carma/hjt_west_1800.txt 1800_west sizedist_1d_plots

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

ntime = nstep/iskip

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
			
			tio2[k,j,i] = float(line[2])/(math.log(ru[j,0]) - math.log(rl[j,0]))
			mass_tio2[k,j,i] = tio2[k,j,i]*(M_TIO2*((4./3.)*np.pi*(r[j,0]*1e-4)**3.))
			
			al2o3[k,j,i] = float(line[3])/(math.log(ru[j,1]) - math.log(rl[j,1]))
			mass_al2o3[k,j,i] = al2o3[k,j,i]*(M_AL2O3*((4./3.)*np.pi*(r[j,1]*1e-4)**3.))
			
			het_fe[k,j,i] = float(line[4])/(math.log(ru[j,2]) - math.log(rl[j,2]))
			mass_het_fe[k,j,i] = het_fe[k,j,i]*(M_FE*((4./3.)*np.pi*(r[j,2]*1e-4)**3.))
			
			het_mg2sio4[k,j,i] = float(line[5])/(math.log(ru[j,3]) - math.log(rl[j,3]))
			mass_het_mg2sio4[k,j,i] = het_mg2sio4[k,j,i]*(M_MG2SIO4*((4./3.)*np.pi*(r[j,3]*1e-4)**3.))
			
			het_cr[k,j,i] = float(line[6])/(math.log(ru[j,4]) - math.log(rl[j,4]))
			mass_het_cr[k,j,i] = cr[k,j,i]*(M_CR*((4./3.)*np.pi*(r[j,4]*1e-4)**3.))
			
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

f_ti = interp2d(r[:,0],p,t_tio2)
f_hfe = interp2d(r[:,2],p,t_het_fe)
f_pfe = interp2d(r[:,7],p,t_pure_fe)
f_mg = interp2d(r[:,3],p,t_het_mg2sio4)
f_hcr = interp2d(r[:,4],p,t_het_cr)
f_pcr = interp2d(r[:,8],p,t_cr)
f_mn = interp2d(r[:,5],p,t_mns)
f_na = interp2d(r[:,6],p,t_na2s)
f_zn = interp2d(r[:,10],p,t_zns)
f_kc = interp2d(r[:,9],p,t_kcl)
f_al = interp2d(r[:,1],p,t_al2o3)

##########################################################################################
##########################################################################################
print 'Plotting it up good!'
##########################################################################################
##########################################################################################
plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['legend.borderpad'] = 0.1
plt.rcParams['legend.labelspacing'] = 0.1
plt.rcParams['legend.handletextpad'] = 0.1
plt.rcParams['font.family'] = 'stixgeneral'
plt.rcParams['font.size'] = 16

radius_array = np.logspace(-4,3,num=100)

h1 = p[k-10]
h2 = p[k-17]
h3 = p[k-20]

plot1 = f_ti(radius_array,h1)+f_hfe(radius_array,h1)+f_pfe(radius_array,h1)+f_mg(radius_array,h1)+f_hcr(radius_array,h1)+f_pcr(radius_array,h1)+f_mn(radius_array,h1)+f_na(radius_array,h1)+f_zn(radius_array,h1)+f_kc(radius_array,h1)+f_al(radius_array,h1) 
plot2 = f_ti(radius_array,h2)+f_hfe(radius_array,h2)+f_pfe(radius_array,h2)+f_mg(radius_array,h2)+f_hcr(radius_array,h2)+f_pcr(radius_array,h2)+f_mn(radius_array,h2)+f_na(radius_array,h2)+f_zn(radius_array,h2)+f_kc(radius_array,h2)+f_al(radius_array,h2) 
plot3 = f_ti(radius_array,h3)+f_hfe(radius_array,h3)+f_pfe(radius_array,h3)+f_mg(radius_array,h3)+f_hcr(radius_array,h3)+f_pcr(radius_array,h3)+f_mn(radius_array,h3)+f_na(radius_array,h3)+f_zn(radius_array,h3)+f_kc(radius_array,h3)+f_al(radius_array,h3)


f,ax = plt.subplots(1,1)

ax.plot(radius_array,plot1, c='b', label=str("{:.1e}".format(h1/1.e6))+' bar', lw=2, alpha = .8)		
ax.plot(radius_array,plot2, c='g', label=str("{:.1e}".format(h2/1.e6))+' bar', lw=2)
ax.plot(radius_array,plot3, c='r', label=str("{:.1e}".format(h3/1.e6))+' bar', lw=2, alpha = .8)	
	

ax.legend(loc='lower left')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(min(radius_array),max(radius_array))
ax.set_ylim(1e-6,3e2)
ax.set_ylabel(r'dN/dLn(r) [cm$^{-3}$]')
ax.set_xlabel(r'Particle Radius [$\mu$m]')


r_normy = np.logspace(-2,2,num=100)
y_normy = lognorm.pdf(r_normy,0.5)

left, bottom, width, height = [0.67, 0.67, 0.2, 0.2]
ax3 = f.add_axes([left, bottom, width, height])
ax3.minorticks_off()
ax3.tick_params(axis='both', which='major', labelsize=12)
ax3.plot(r_normy,y_normy*1e1, c='k')
ax3.set_xlim(1e-4,1e4)
ax3.set_ylim(1e-6,1e2)
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xticks([1e-4, 1e-2, 1e0,1e2,1e4])
ax3.set_yticks([1e-6,1e-4,1e-2,1e0,1e2])
ax3.tick_params(axis='both',which='minor',bottom='off', top ='off' , left ='off' , right ='off' )


plt.savefig(path+'/nd_dist_'+name+'.jpg', format = 'jpg', bbox_inches='tight')
plt.close()


