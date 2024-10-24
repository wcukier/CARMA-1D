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

'''
print r[:,0], 'TiO2'
print r[:,1], 'Fe'
print r[:,2], 'Mg'
print r[:,5], "Cr"
print r[:,6], 'Mn'
print r[:,7], 'Na'
print r[:,8], 'Zn'
print r[:,9], 'K'
print r[:,10], 'Al'
'''
			
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

data.close()


f_tio2 = interp2d(r[:,0],p/1.e6,t_tio2)
f_het_fe = interp2d(r[:,1],p/1.e6,t_het_fe) 
f_het_mg2sio4 = interp2d(r[:,1],p/1.e6,t_het_mg2sio4)
f_pure_fe = interp2d(r[:,0],p/1.e6,t_pure_fe)
f_pure_mg2sio4 = interp2d(r[:,0],p/1.e6,t_pure_mg2sio4)
f_cr = interp2d(r[:,0],p/1.e6,t_cr)
f_mns = interp2d(r[:,0],p/1.e6,t_mns)
f_na2s = interp2d(r[:,0],p/1.e6,t_na2s)
f_zns = interp2d(r[:,0],p/1.e6,t_zns)
f_kcl = interp2d(r[:,0],p/1.e6,t_kcl)
f_al2o3 = interp2d(r[:,0],p/1.e6,t_al2o3)

wavelength,extinct_ti, ssa_ti, scat_ti = np.loadtxt('rf_data/ti_rad_1.txt', unpack=True)

extinct_tio2 = zeros((len(wavelength),nbin))
extinct_fe = zeros((len(wavelength),nbin))
extinct_mg2sio4 = zeros((len(wavelength),nbin))
extinct_cr = zeros((len(wavelength),nbin))
extinct_mns = zeros((len(wavelength),nbin))
extinct_na2s = zeros((len(wavelength),nbin))
extinct_zns = zeros((len(wavelength),nbin))
extinct_kcl = zeros((len(wavelength),nbin))
extinct_al2o3 = zeros((len(wavelength),nbin))

assym_tio2 = zeros((len(wavelength),nbin))
assym_fe = zeros((len(wavelength),nbin))
assym_mg2sio4 = zeros((len(wavelength),nbin))
assym_cr = zeros((len(wavelength),nbin))
assym_mns = zeros((len(wavelength),nbin))
assym_na2s = zeros((len(wavelength),nbin))
assym_zns = zeros((len(wavelength),nbin))
assym_kcl = zeros((len(wavelength),nbin))
assym_al2o3 = zeros((len(wavelength),nbin))

scat_tio2 = zeros((len(wavelength),nbin))
scat_fe = zeros((len(wavelength),nbin))
scat_mg2sio4 = zeros((len(wavelength),nbin))
scat_cr = zeros((len(wavelength),nbin))
scat_mns = zeros((len(wavelength),nbin))
scat_na2s = zeros((len(wavelength),nbin))
scat_zns = zeros((len(wavelength),nbin))
scat_kcl = zeros((len(wavelength),nbin))
scat_al2o3 = zeros((len(wavelength),nbin))

ssa_ti = zeros((len(wavelength),nbin))
ssa_fe = zeros((len(wavelength),nbin))
ssa_mg = zeros((len(wavelength),nbin))
ssa_cr = zeros((len(wavelength),nbin))
ssa_mn = zeros((len(wavelength),nbin))
ssa_na = zeros((len(wavelength),nbin))
ssa_zn = zeros((len(wavelength),nbin))
ssa_k = zeros((len(wavelength),nbin))
ssa_al = zeros((len(wavelength),nbin))

for j in range(nbin):
	wi, exti, assymti, scatti = np.loadtxt('rf_data/ti_rad_'+str(j)+'.txt', unpack=True)
	wi, exmg, assymmg, scatmg = np.loadtxt('rf_data/mg_rad_'+str(j)+'.txt', unpack=True)
	wi, exfe, assymfe, scatfe = np.loadtxt('rf_data/fe_rad_'+str(j)+'.txt', unpack=True)
	wi, excr, assymcr, scatcr = np.loadtxt('rf_data/cr_rad_'+str(j)+'.txt', unpack=True)
	wi, exmn, assymmn, scatmn = np.loadtxt('rf_data/mn_rad_'+str(j)+'.txt', unpack=True)
	wi, exna, assymna, scatna = np.loadtxt('rf_data/na_rad_'+str(j)+'.txt', unpack=True)
	wi, exzn, assymzn, scatzn = np.loadtxt('rf_data/zn_rad_'+str(j)+'.txt', unpack=True)
	wi, exk, assymk, scatk = np.loadtxt('rf_data/k_rad_'+str(j)+'.txt', unpack=True)
	wi, exal, assymal, scatal = np.loadtxt('rf_data/al_rad_'+str(j)+'.txt', unpack=True)
	
	extinct_tio2[:,j] = exti
	extinct_mg2sio4[:,j] = exmg
	extinct_fe[:,j] = exfe
	extinct_cr[:,j] = excr
	extinct_mns[:,j] = exmn
	extinct_na2s[:,j] = exna
	extinct_zns[:,j] = exzn
	extinct_kcl[:,j] = exk
	extinct_al2o3[:,j] = exal
	
	assym_tio2[:,j] = assymti
	assym_mg2sio4[:,j] = assymmg
	assym_fe[:,j] = assymfe
	assym_cr[:,j] = assymcr
	assym_mns[:,j] = assymmn
	assym_na2s[:,j] = assymna
	assym_zns[:,j] = assymzn
	assym_kcl[:,j] = assymk
	assym_al2o3[:,j] = assymal
	
	scat_tio2[:,j] = scatti
	scat_mg2sio4[:,j] = scatmg
	scat_fe[:,j] = scatfe
	scat_cr[:,j] = scatcr
	scat_mns[:,j] = scatmn
	scat_na2s[:,j] = scatna
	scat_zns[:,j] = scatzn
	scat_kcl[:,j] = scatk
	scat_al2o3[:,j] = scatal
	
	ssa_ti[:,j] = scatti/exti
	ssa_mg[:,j] = scatmg/exmg
	ssa_fe[:,j] = scatfe/exfe
	ssa_cr[:,j] = scatcr/excr
	ssa_mn[:,j] = scatmn/exmn
	ssa_na[:,j] = scatna/exna
	ssa_zn[:,j] = scatzn/exzn
	ssa_k[:,j] = scatk/exk
	ssa_al[:,j] = scatal/exal
	
	ssa_ti[:,j][ssa_ti[:,j]> 1.] = 1.
	ssa_mg[:,j][ssa_mg[:,j]> 1.] = 1.
	ssa_fe[:,j][ssa_fe[:,j]> 1.] = 1.
	ssa_cr[:,j][ssa_cr[:,j]> 1.] = 1.
	ssa_mn[:,j][ssa_mn[:,j]> 1.] = 1.
	ssa_na[:,j][ssa_na[:,j]> 1.] = 1.
	ssa_zn[:,j][ssa_zn[:,j]> 1.] = 1.
	ssa_k[:,j][ssa_k[:,j]> 1.] = 1.
	ssa_al[:,j][ssa_al[:,j]> 1.] = 1.


print ssa_ti[:,:], 'ssa'

grid_data = np.loadtxt('sonora_profiles/t'+name+'g100f1_m0.0.dat', skiprows = 1)
p_grid = grid_data.T[1]

f_height = interp1d(p/1.e6,z)

total_opacity_int = np.zeros((len(p_grid), len(wavelength), nbin))
total_scattering_int = np.zeros((len(p_grid), len(wavelength), nbin))
total_assym_int = np.zeros((len(p_grid), len(wavelength), nbin))
SSA_int = np.zeros((len(p_grid), len(wavelength), nbin))

total_opacity = np.zeros((len(p_grid), len(wavelength)))
total_scattering = np.zeros((len(p_grid), len(wavelength)))
total_assym = np.zeros((len(p_grid), len(wavelength)))

SSA = np.zeros((len(p_grid), len(wavelength)))

for idx, pp in enumerate(p_grid):
	next = int(idx+1)
	for iw, w in enumerate(wavelength):
		for j in range(nbin):
			if next > 90 or pp > max(p)/1.e6 or p_grid[next] > max(p)/1.e6 or idx ==0:
				total_opacity_int[idx, iw, j] = 0.
				total_scattering_int[idx, iw, j] = 0.
				
				frac_ti = 0.
				frac_mg = 0.
				frac_fe = 0.
				frac_cr = 0.
				frac_mn = 0.
				frac_na = 0.
				frac_al = 0.
				frac_zn = 0.
				frac_k = 0.
								
			elif idx == 1:
				id = int(90)
				dz = (max(z)-f_height(p_grid[next]))*100.
				
				opac_ti = t_tio2[id,j]*extinct_tio2[iw,j]
				opac_mg = t_het_mg2sio4[id,j]*extinct_mg2sio4[iw,j]
				opac_fe = t_het_fe[id,j]*extinct_fe[iw,j]
				opac_cr = t_cr[id,j]*extinct_cr[iw,j]
				opac_mn = t_mns[id,j]*extinct_mns[iw,j]
				opac_zn = t_zns[id,j]*extinct_zns[iw,j]
				opac_na = t_na2s[id,j]*extinct_na2s[iw,j]
				opac_al = t_al2o3[id,j]*extinct_al2o3[iw,j]
				opac_k = t_kcl[id,j]*extinct_kcl[iw,j]
				
				sca_ti = t_tio2[id,j]*scat_tio2[iw,j]
				sca_mg = t_het_mg2sio4[id,j]*scat_mg2sio4[iw,j]
				sca_fe = t_het_fe[id,j]*scat_fe[iw,j]
				sca_cr = t_cr[id,j]*scat_cr[iw,j]
				sca_mn = t_mns[id,j]*scat_mns[iw,j]
				sca_zn = t_zns[id,j]*scat_zns[iw,j]
				sca_na = t_na2s[id,j]*scat_na2s[iw,j]
				sca_al = t_al2o3[id,j]*scat_al2o3[iw,j]
				sca_k = t_kcl[id,j]*scat_kcl[iw,j]
				
				total_opacity_int[idx, iw, j] = (opac_ti+opac_mg+opac_fe+opac_cr+opac_mn+opac_zn+opac_na+opac_al+opac_k)*dz
				total_scattering_int[idx, iw, j] = (sca_ti+sca_mg+sca_fe+sca_cr+sca_mn+sca_zn+sca_na+sca_al+sca_k)*dz
				
				frac_ti = opac_ti/(opac_ti+opac_mg+opac_fe+opac_cr+opac_mn+opac_zn+opac_na+opac_al+opac_k)
				frac_mg = opac_mg/(opac_ti+opac_mg+opac_fe+opac_cr+opac_mn+opac_zn+opac_na+opac_al+opac_k)
				frac_fe = opac_fe/(opac_ti+opac_mg+opac_fe+opac_cr+opac_mn+opac_zn+opac_na+opac_al+opac_k)
				frac_cr = opac_cr/(opac_ti+opac_mg+opac_fe+opac_cr+opac_mn+opac_zn+opac_na+opac_al+opac_k)
				frac_mn = opac_mn/(opac_ti+opac_mg+opac_fe+opac_cr+opac_mn+opac_zn+opac_na+opac_al+opac_k)
				frac_na = opac_na/(opac_ti+opac_mg+opac_fe+opac_cr+opac_mn+opac_zn+opac_na+opac_al+opac_k)
				frac_al = opac_al/(opac_ti+opac_mg+opac_fe+opac_cr+opac_mn+opac_zn+opac_na+opac_al+opac_k)
				frac_zn = opac_zn/(opac_ti+opac_mg+opac_fe+opac_cr+opac_mn+opac_zn+opac_na+opac_al+opac_k)
				frac_k = opac_k/(opac_ti+opac_mg+opac_fe+opac_cr+opac_mn+opac_zn+opac_na+opac_al+opac_k)
				
			else:
				dz = (f_height(pp)-f_height(p_grid[next]))*100.
				
				opac_ti = f_tio2(r[j,1], pp)*extinct_tio2[iw,j]
				opac_mg = f_het_mg2sio4(r[j,1], pp)*extinct_mg2sio4[iw,j]
				opac_fe = f_het_fe(r[j,1], pp)*extinct_fe[iw,j]
				opac_cr = f_cr(r[j,1], pp)*extinct_cr[iw,j]
				opac_mn = f_mns(r[j,1], pp)*extinct_mns[iw,j]
				opac_zn = f_zns(r[j,1], pp)*extinct_zns[iw,j]
				opac_na = f_na2s(r[j,1], pp)*extinct_na2s[iw,j]
				opac_al = f_al2o3(r[j,1], pp)*extinct_al2o3[iw,j]
				opac_k = f_kcl(r[j,1], pp)*extinct_kcl[iw,j]
				
				sca_ti = f_tio2(r[j,1], pp)*scat_tio2[iw,j]
				sca_mg = f_het_mg2sio4(r[j,1], pp)*scat_mg2sio4[iw,j]
				sca_fe = f_het_fe(r[j,1], pp)*scat_fe[iw,j]
				sca_cr = f_cr(r[j,1], pp)*scat_cr[iw,j]
				sca_mn = f_mns(r[j,1], pp)*scat_mns[iw,j]
				sca_zn = f_zns(r[j,1], pp)*scat_zns[iw,j]
				sca_na = f_na2s(r[j,1], pp)*scat_na2s[iw,j]
				sca_al = f_al2o3(r[j,1], pp)*scat_al2o3[iw,j]
				sca_k = f_kcl(r[j,1], pp)*scat_kcl[iw,j]
				
				total_opacity_int[idx, iw, j] = (opac_ti+opac_mg+opac_fe+opac_cr+opac_mn+opac_zn+opac_na+opac_al+opac_k)*dz
				total_scattering_int[idx, iw, j] = (sca_ti+sca_mg+sca_fe+sca_cr+sca_mn+sca_zn+sca_na+sca_al+sca_k)*dz
									   
				frac_ti = opac_ti/(opac_ti+opac_mg+opac_fe+opac_cr+opac_mn+opac_zn+opac_na+opac_al+opac_k)
				frac_mg = opac_mg/(opac_ti+opac_mg+opac_fe+opac_cr+opac_mn+opac_zn+opac_na+opac_al+opac_k)
				frac_fe = opac_fe/(opac_ti+opac_mg+opac_fe+opac_cr+opac_mn+opac_zn+opac_na+opac_al+opac_k)
				frac_cr = opac_cr/(opac_ti+opac_mg+opac_fe+opac_cr+opac_mn+opac_zn+opac_na+opac_al+opac_k)
				frac_mn = opac_mn/(opac_ti+opac_mg+opac_fe+opac_cr+opac_mn+opac_zn+opac_na+opac_al+opac_k)
				frac_na = opac_na/(opac_ti+opac_mg+opac_fe+opac_cr+opac_mn+opac_zn+opac_na+opac_al+opac_k)
				frac_al = opac_al/(opac_ti+opac_mg+opac_fe+opac_cr+opac_mn+opac_zn+opac_na+opac_al+opac_k)
				frac_zn = opac_zn/(opac_ti+opac_mg+opac_fe+opac_cr+opac_mn+opac_zn+opac_na+opac_al+opac_k)
				frac_k  =  opac_k/(opac_ti+opac_mg+opac_fe+opac_cr+opac_mn+opac_zn+opac_na+opac_al+opac_k)
			
			total_assym_int[idx,iw,j] = (assym_mg2sio4[iw,j]*frac_mg+assym_tio2[iw,j]*frac_ti+assym_fe[iw,j]*frac_fe+assym_cr[iw,j]*frac_cr+assym_mns[iw,j]*frac_mn+
										assym_na2s[iw,j]*frac_na+assym_al2o3[iw,j]*frac_al+assym_zns[iw,j]*frac_zn+assym_kcl[iw,j]*frac_k)
			
			SSA_int[idx,iw,j] = (ssa_mg[iw,j]*frac_mg + ssa_ti[iw,j]*frac_ti + ssa_fe[iw,j]*frac_fe + ssa_cr[iw,j]*frac_cr + ssa_mn[iw,j]*frac_mn +
								ssa_na[iw,j]*frac_na + ssa_al[iw,j]*frac_al + ssa_zn[iw,j]*frac_zn + ssa_k[iw,j]*frac_k)
			
									   
		total_opacity[idx,iw] = sum(total_opacity_int[idx,iw,:])
		total_scattering[idx,iw] = sum(total_scattering_int[idx,iw,:])
		
		for j in range(nbin):
			frac_bin	= total_opacity_int[idx,iw,j]/total_opacity[idx,iw]
			total_assym_int[idx,iw,j] = total_assym_int[idx,iw,j]*frac_bin
			SSA_int[idx,iw,j] = SSA_int[idx,iw,j]*frac_bin
		
		total_assym[idx,iw] = sum(total_assym_int[idx,iw,:])	
		SSA[idx,iw] = sum(SSA_int[idx,iw,:])
		
		if math.isnan(SSA[idx,iw]):
			SSA[idx,iw] = 0.
		if math.isnan(total_assym[idx,iw]):
			total_assym[idx,iw] = 0.
		if math.isnan(total_opacity[idx,iw]):
			total_opacity[idx,iw] = 0.
		
		total_opacity[idx,iw] = np.abs(total_opacity[idx,iw])
		total_assym[idx,iw] = np.abs(total_assym[idx,iw])
		SSA[idx,iw] = np.abs(SSA[idx,iw])
		
	total_assym[idx,:][total_assym[idx,:]> 1.] = 1.
	SSA[idx,:][SSA[idx,:]> 1.] = 1.
	

for idx, pp in enumerate(p_grid):
	print max(SSA[idx,:]), 'SSA'
	print max(total_opacity[idx,:]), 'Opaque'
	print max(total_assym[idx,:]), 'Assym' 

def to_precision(x,p):
    """
    returns a string representation of x formatted with a precision of p

    Based on the webkit javascript implementation taken from here:
    https://code.google.com/p/webkit-mirror/source/browse/JavaScriptCore/kjs/number_object.cpp
    """

    x = float(x)

    if x == 0.:
        return "0." + "0"*(p-1)

    out = []

    if x < 0:
        out.append("-")
        x = -x

    e = int(math.log10(x))
    tens = math.pow(10, e - p + 1)
    n = math.floor(x/tens)

    if n < math.pow(10, p - 1):
        e = e -1
        tens = math.pow(10, e - p+1)
        n = math.floor(x / tens)

    if abs((n + 1.) * tens - x) <= abs(n * tens -x):
        n = n + 1

    if n >= math.pow(10,p):
        n = n / 10.
        e = e + 1

    m = "%.*g" % (p, n)

    if e < -2 or e >= p:
        out.append(m[0])
        if p > 1:
            out.append(".")
            out.extend(m[1:p])
        out.append('e')
        if e > 0:
            out.append("+")
        out.append(str(e))
    elif e == (p -1):
        out.append(m)
    elif e >= 0:
        out.append(m[:e+1])
        if e+1 < len(m):
            out.append(".")
            out.extend(m[e+1:])
    else:
        out.append("0.")
        out.extend(["0"]*-(e+1))
        out.append(m)

    return "".join(out)

column1 = []
column2 = []
column3 = []
column4 = []
column5 = []

for idx, pp in enumerate(p_grid):
	for iw, w in enumerate(wavelength):
		column1.append(int(idx+1))
		column2.append(int(iw+1))
		column3.append(total_opacity[idx,iw])
		column4.append(to_precision(total_assym[idx,iw],8))
		column5.append(to_precision(SSA[idx,iw],8))

np.savetxt('t'+name+'g100f1_m0.0_CARMAout.cld', np.transpose([column1, column2, column3, column4, column5]), fmt='%s')

'''
outfile = open('t'+name+'g100f1_m0.0_CARMAout.cld', 'w')
for idx, pp in enumerate(p_grid):
	for iw, w in enumerate(wavelength):
		outfile.write('%.2f  %.4f\n' % (int(idx+1), int(iw+1), total_opacity[idx,iw], o_precision(total_assym[idx,iw],8), to_precision(SSA[idx,iw],8)))
		
'''




















