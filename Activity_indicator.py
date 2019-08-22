"""
Created on 16 Aug 2019

@author: jzhao
"""

import sys
import shutil
from astropy.io import fits
import glob
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
from scipy import stats



#############################################################
###################### Define function ######################
#############################################################

def read_rdb(filename):
    
    f = open(filename, 'r')
    data = f.readlines()
    f.close()
    
    z=0
    while data[z][:2] == '# ' or data[z][:2] == ' #':
        z += 1

    key = str.split(data[z+0][:-1],'\t')
    output = {}
    for i in range(len(key)): output[key[i]] = []
    
    for line in data[z+2:]:
        qq = str.split(line[:-1],'\t')
        for i in range(len(key)):
            try: value = float(qq[i])
            except ValueError: value = qq[i]
            output[key[i]].append(value)

    return output


def gaussian(x, a, mu, sigma, C):
    val = a / ( (2*math.pi)**0.5 * sigma ) * np.exp( -(x - mu)**2 / (2*sigma**2) ) + C
    return val

#################################
# Results from the HARPS header #
#################################
data_ccf    = read_rdb('/Volumes/DataSSD/SOAP_2/outputs/02.01/phase_rv.dat')
RV_HARPS 	= np.array(data_ccf['RV_tot']) * 1000   	# m/s
FWHM_HARPS 	= np.array(data_ccf['Fwhm_tot'])   			# km/s
BIS_HARPS 	= np.array(data_ccf['Bis_Span_tot'])*1000   # m/s
BIS_HARPS 	= BIS_HARPS - np.mean(BIS_HARPS)

#############################################
# Define variables 
# NOISE 	= 5000
# NOISE 	= 2000
NOISE 	= False

if NOISE:
	N = 10000			# number of files 
else:
	N = 100
RV 		= np.zeros(N)
FWHM	= np.zeros(N)
V_span	= np.zeros(N)
BIS 	= np.zeros(N)


v_CCF 	= np.linspace(-20, 20, 401)
idx0 	= (v_CCF >= -10) & (v_CCF <= 10)
v_ccf  	= v_CCF[idx0]

for i in range(N):
	if NOISE:
		filename	= '/Volumes/DataSSD/SOAP_2/outputs/02.01/CCF_dat_SN=' + str(NOISE) + '/CCF' + str(i) + '.dat'
	else:
		filename 	= '/Volumes/DataSSD/SOAP_2/outputs/02.01/CCF_dat/CCF' + str(i) + '.dat'
	CCF 		= 1 - np.loadtxt(filename)
	ccf 		= CCF[idx0]
	# popt, pcov  = curve_fit(gaussian, v_ccf, ccf)
	popt, pcov  = curve_fit(gaussian, v_CCF, CCF)
	popt0 		= popt 								# Keep a snapshot
	RV[i]   	= popt[1]    						# km/s
	FWHM[i]		= popt[2] * math.sqrt(8*math.log(2))
	sigma 		= FWHM[i]/2

	###################
	# Calculate Vspan #
	###################
	# σ is the width of the CCF
	# I take σ as the FWHM

	# the upper part of the CCF is defined in the range [–∞:–1σ][+1σ:+∞] 
	idx     = (v_CCF < (RV[i] - sigma)) + (v_CCF > (RV[i] + sigma))
	x_span  = v_CCF[idx]
	y_span  = CCF[idx]/CCF.max()
	popt, pcov  = curve_fit(gaussian, x_span, y_span, popt0)
	V_high = popt[1]

	# the lower part is defined in the range given by [–∞:-3σ][–1σ:+1σ][+3σ:+∞].
	idx     = (v_CCF < (RV[i] - 3*sigma)) + (v_CCF > (RV[i] + 3*sigma)) + (v_CCF > (RV[i] - sigma)) * (v_CCF < (RV[i] + sigma))
	x_span  = v_CCF[idx]
	y_span  = CCF[idx]/CCF.max()
	popt, pcov  = curve_fit( gaussian, x_span, y_span, popt0)
	V_low = popt[1]
	V_span[i] = (V_high - V_low) * 1000		# in unit of m/s

	###########################
	# Calculate bisector span #
	###########################
	ccf 	= ccf / (popt0[0] / ( (2*math.pi)**0.5 * popt0[2] )) 	# normalization
	idx_l 	= v_ccf < 0
	idx_r 	= v_ccf > 0
	v_l 	= v_ccf[idx_l]
	v_r 	= v_ccf[idx_r]
	ccf_l 	= ccf[idx_l]
	ccf_r 	= ccf[idx_r]
	idx_ll 	= (ccf_l <= 0.9) & (ccf_l >= 0.6)
	idx_rr 	= (ccf_r <= 0.9) & (ccf_r >= 0.6)

	# rotate CCF for interpolation
	x 		= np.linspace(0.6, 0.9, 101)
	cs 		= CubicSpline(ccf_r[idx_rr][::-1], v_r[idx_rr][::-1])
	temp1 	= np.mean(cs(x)) 		# mean of v_r at the bottom
	cs 		= CubicSpline(ccf_l[idx_ll], v_l[idx_ll])
	temp2 	= np.mean(cs(x)) 		# mean of v_l at the bottom
	v_bottom= (temp1 + temp2)/2  	# Note that the line profile is reversed

	# plt.plot(v_r[idx_rr], ccf_r[idx_rr], '.'); plt.show()
	# plt.plot(ccf_r[idx_rr], v_r[idx_rr], '.'); plt.show()

	idx_ll 	= (ccf_l <= 0.4) & (ccf_l >= 0.1)
	idx_rr 	= (ccf_r <= 0.4) & (ccf_r >= 0.1)
	x 		= np.linspace(0.1, 0.4, 101)
	cs 		= CubicSpline(ccf_r[idx_rr][::-1], v_r[idx_rr][::-1])
	temp1 	= np.mean(cs(x)) 		# mean of v_r at the bottom
	cs 		= CubicSpline(ccf_l[idx_ll], v_l[idx_ll])
	temp2 	= np.mean(cs(x)) 		# mean of v_l at the bottom
	v_top 	= (temp1 + temp2)/2  	# Note that the line profile is reversed

	BIS[i] 	= v_top - v_bottom

BIS = (BIS - np.mean(BIS))*1000	# m/s


#####################
# Correlation plots #
#####################
RV 		= (RV - np.mean(RV))*1000
V_span 	= V_span - np.mean(V_span)
if NOISE:
	DIR 		= '/Volumes/DataSSD/MATLAB_codes/Project180316-shfit_in_FT/jitter_SN=' + str(NOISE) +'/'
	SAVE_NAME 	= DIR + 'correlation_SN' + str(NOISE) + '.png'
	alpha = 0.01
else:
	DIR 		= '/Volumes/DataSSD/MATLAB_codes/Project180316-shfit_in_FT/jitter_noise_free/'
	SAVE_NAME 	= DIR + 'correlation_noise_free.png'
	alpha = 0.5
dRV_L 	= np.loadtxt(DIR + 'YY.txt')
dRV_H 	= np.loadtxt(DIR + 'ZZ.txt')

left  = 0.08  # the left side of the subplots of the figure
right = 0.95    # the right side of the subplots of the figure
bottom = 0.2   # the bottom of the subplots of the figure
top = 0.8      # the top of the subplots of the figure
wspace = 0.6   # the amount of width reserved for blank space between subplots
hspace = 0.2   # the amount of height reserved for white space between subplots
Nx = 3
Ny = 3

plt.rcParams.update({'font.size': 20})
fig, axes = plt.subplots(figsize=(20, 4))
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

axes_1 = plt.subplot(151)
plt.plot(RV, FWHM, 'ko', alpha=alpha)
plt.xlabel('Jitter [m/s]')
plt.ylabel('FWHM [km/s]')
r1, p = stats.pearsonr(RV, FWHM)
plt.title(r'$\rho = {0:.2f}$'.format(r1))
axes_1.xaxis.set_major_locator(plt.MaxNLocator(Nx))
axes_1.yaxis.set_major_locator(plt.MaxNLocator(Ny))

axes_2 = plt.subplot(152)
plt.plot(RV, BIS, 'ko', alpha=alpha)
plt.xlabel('Jitter [m/s]')
plt.ylabel('BIS [m/s]')
r1, p = stats.pearsonr(RV, BIS)
plt.title(r'$\rho = {0:.2f}$'.format(r1))
axes_2.xaxis.set_major_locator(plt.MaxNLocator(Nx))
axes_2.yaxis.set_major_locator(plt.MaxNLocator(Ny))

axes_3 = plt.subplot(153)
plt.plot(RV, V_span, 'ko', alpha=alpha)
plt.xlabel('Jitter [m/s]')
plt.ylabel(r'$V_{span}$ [m/s]')
r1, p = stats.pearsonr(RV, V_span)
plt.title(r'$\rho = {0:.2f}$'.format(r1))
axes_3.xaxis.set_major_locator(plt.MaxNLocator(Nx))
axes_3.yaxis.set_major_locator(plt.MaxNLocator(Ny))

axes_4 = plt.subplot(154)
plt.plot(RV, dRV_L, 'ko', alpha=alpha)
plt.xlabel('Jitter [m/s]')
plt.ylabel(r'$\Delta RV_L$ [m/s]')
r1, p = stats.pearsonr(RV, dRV_L)
plt.title(r'$\rho = {0:.2f}$'.format(r1))
axes_4.xaxis.set_major_locator(plt.MaxNLocator(Nx))
axes_4.yaxis.set_major_locator(plt.MaxNLocator(Ny))

axes_5 = plt.subplot(155)
plt.plot(RV, dRV_H, 'ko', alpha=alpha)
plt.xlabel('Jitter [m/s]')
plt.ylabel(r'$\Delta RV_H$ [m/s]')
r1, p = stats.pearsonr(RV, dRV_H)
plt.title(r'$\rho = {0:.2f}$'.format(r1))
axes_5.xaxis.set_major_locator(plt.MaxNLocator(Nx))
axes_5.yaxis.set_major_locator(plt.MaxNLocator(Ny))

plt.savefig(SAVE_NAME)
plt.show()


######################################################
# Testing -- how consistent with HARPS measurements? #
######################################################
if 0:
	# BIS #
	fig = plt.figure()
	frame1 = fig.add_axes((0.15, 0.32, 0.8, 0.6)) # left, bottom, width, height
	plt.plot(BIS_HARPS, BIS, '.')
	plt.ylabel("my BIS [m/s]")	
	plt.ylim(-1.5, 1.5)
	frame1.set_xticklabels([])
	plt.title('BIS')

	fig.add_axes((.15, .1, .8, .2))  # left, bottom, width, height
	plt.xlabel(r"$BIS_{HARPS}$ [m/s]")	
	plt.ylabel("Residual [m/s]")	
	plt.ylim(-0.5, 0.5)
	plt.plot(BIS_HARPS, BIS-BIS_HARPS, '.')
	plt.show()

	# FWHM #
	fig = plt.figure()
	frame1 = fig.add_axes((0.15, 0.32, 0.8, 0.6)) # left, bottom, width, height
	plt.plot(FWHM_HARPS, FWHM, '.')
	plt.ylabel("my FWHM [km/s]")	
	# plt.ylim(6.4415, 6.4475)
	frame1.set_xticklabels([])
	plt.title('FWHM')

	fig.add_axes((.15, .1, .8, .2))  # left, bottom, width, height
	plt.xlabel(r"$FWHM_{HARPS}$ [km/s]")	
	plt.ylabel("Residual [km/s]")	
	plt.ylim(-0.001, 0.001)
	plt.plot(FWHM_HARPS, FWHM-FWHM_HARPS, '.')
	plt.show()

if 0:
	plt.plot(v_CCF, CCF); plt.show()
	# best fits 
	plt.plot(v_ccf, ccf) 
	# plt.plot(v_ccf, gaussian(v_ccf, *popt))
	plt.show()

	# RV, FWHM
	plt.plot(RV_HARPS, RV-RV_HARPS); plt.show()
	plt.plot(FWHM_HARPS, FWHM-FWHM_HARPS); plt.show()
	plt.plot(FWHM_HARPS, V_span); plt.show()
	
	# V_span
	plt.plot(v_CCF, gaussian(v_CCF, *popt))
	plt.plot(x_span, y_span, '.') 
	plt.show()


