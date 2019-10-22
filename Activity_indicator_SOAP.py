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


###################
# spot/plage size #
###################
N_sampling  = 20
S1          = np.linspace(0.5, 10, num=N_sampling, endpoint=True)/100 # size in area
r1 		= np.zeros((5, N_sampling))
k 		= np.zeros((5, N_sampling))

for n in range(N_sampling):
	size = np.sqrt(2*S1[n])
	folder_name = 'CCF_PROT=25.05_i=90.00_lon=(180.0,0.0,0.0,0.0)_lat=(0.0,0.0,0.0,0.0)_size=(' + '%.4f' % round(size, 4) + ',0.0000,0.0000,0.0000)'
	DIR     = '/Volumes/DataSSD/SOAP_2/outputs/' + folder_name

##################
# spot/plage lat #
##################
# N_sampling  = 18
# lat         = np.linspace(0, 85, num=N_sampling, endpoint=True)
# r1 		= np.zeros((5, N_sampling))
# k 		= np.zeros((5, N_sampling))

# for n in range(N_sampling):
# 	folder_name = 'CCF_PROT=25.05_i=90.00_lon=(180.0,0.0,0.0,0.0)_lat=(' + '%.1f' % lat[n] +  ',0.0,0.0,0.0)_size=(0.1000,0.0000,0.0000,0.0000)'
# 	DIR     = '/Volumes/DataSSD/SOAP_2/outputs/' + folder_name

######
# SN #
######
# N_sampling  = 10
# SN         	= np.linspace(1000, 10000, num=N_sampling, endpoint=True)
# r1 			= np.zeros((5, N_sampling))

# for n in range(N_sampling):
# 	folder_name = 'CCF_PROT=25.05_i=90.00_lon=(180.0,0.0,0.0,0.0)_lat=(0.0,0.0,0.0,0.0)_size=(0.1000,0.0000,0.0000,0.0000)'
# 	DIR     = '/Volumes/DataSSD/SOAP_2/outputs/' + folder_name

	data_ccf    = read_rdb(DIR + '/' + folder_name + '.rdb')
	RV_HARPS 	= np.array(data_ccf['RV_tot']) * 1000   	# m/s
	FWHM_HARPS 	= np.array(data_ccf['Fwhm_tot'])   			# km/s
	BIS_HARPS 	= np.array(data_ccf['Bis_Span_tot'])*1000   # m/s
	BIS_HARPS 	= BIS_HARPS - np.mean(BIS_HARPS)

	#############################################

	N = 100
	RV 		= np.zeros(N)
	FWHM	= np.zeros(N)
	V_span	= np.zeros(N)
	BIS 	= np.zeros(N)


	v_CCF 	= np.linspace(-20, 20, 401)
	idx0 	= (v_CCF >= -10) & (v_CCF <= 10)
	v_ccf  	= v_CCF[idx0]

	for i in range(N):
		filename	= DIR + '/fits/CCF' + str(i) + '.dat'
		# filename	= DIR + '/CCF_noise' + '%d' %SN[n] + '/' + str(i) + '.dat'
		CCF 		= 1 - np.loadtxt(filename)
		ccf 		= CCF[idx0]
		# popt, pcov  = curve_fit(gaussian, v_ccf, ccf)
		popt, pcov  = curve_fit(gaussian, v_CCF, CCF)
		popt0 		= popt 								# Keep a snapshot
		RV[i]   	= popt[1]    						# km/s
		FWHM[i]		= popt[2] * math.sqrt(8*math.log(2))
		sigma 		= FWHM[i]/2 						# performs better than sigma = popt[2]
		# sigma 		= popt[2]

		###################
		# Calculate Vspan #
		###################
		# σ is the width of the CCF
		# I take σ as the FWHM

		# oversample v_CCF
		from scipy.interpolate import interp1d
		f 		= interp1d(v_CCF, CCF, kind='cubic')
		v_CCFi 	= np.linspace(-10, 10, 2000)
		CCFi 	= f(v_CCFi)

		# the upper part of the CCF is defined in the range [–∞:–1σ][+1σ:+∞] 
		idx     = (v_CCFi < (RV[i] - sigma)) + (v_CCFi > (RV[i] + sigma))
		x_span  = v_CCFi[idx]
		y_span  = CCFi[idx]/CCFi.max()
		popt, pcov  = curve_fit(gaussian, x_span, y_span, popt0)
		V_high = popt[1]

		# the lower part is defined in the range given by [–∞:-3σ][–1σ:+1σ][+3σ:+∞].
		idx     = (v_CCFi < (RV[i] - 3*sigma)) + (v_CCFi > (RV[i] + 3*sigma)) + (v_CCFi > (RV[i] - sigma)) * (v_CCFi < (RV[i] + sigma))
		x_span  = v_CCFi[idx]
		y_span  = CCFi[idx]/CCFi.max()
		popt, pcov  = curve_fit( gaussian, x_span, y_span, popt0)
		# plt.plot(x_span, y_span, '.')
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
	# plt.show()

	#####################
	# Correlation plots #
	#####################
	RV 		= (RV - np.mean(RV))*1000
	V_span 	= V_span - np.mean(V_span)
	SAVE_NAME 	= DIR + '/correlation' + '%.4f' % round(size, 4) + '.png'
	# SAVE_NAME 	= DIR + '/correlation' + '%.1f' % lat[n] + '.png'
	# SAVE_NAME 		= DIR + '/correlation' + str(SN[n]) + '.png'
	alpha = 0.5
	dRV_L 	= np.loadtxt(DIR + '/YY.txt')
	dRV_H 	= np.loadtxt(DIR + '/ZZ.txt')
	# dRV_L 	= np.loadtxt(DIR + '/CCF_noise' + str(SN[n]) + '/' + '/YY.txt')
	# dRV_H 	= np.loadtxt(DIR + '/CCF_noise' + str(SN[n]) + '/' + '/ZZ.txt')	

	left  	= 0.08  # the left side of the subplots of the figure
	right 	= 0.95    # the right side of the subplots of the figure
	bottom 	= 0.2   # the bottom of the subplots of the figure
	top 	= 0.8      # the top of the subplots of the figure
	wspace 	= 0.6   # the amount of width reserved for blank space between subplots
	hspace 	= 0.2   # the amount of height reserved for white space between subplots
	fontsize=18
	w 		= np.ones(N)
	Nx 		= 3
	Ny 		= 3

	plt.rcParams.update({'font.size': 20})
	fig, axes = plt.subplots(figsize=(20, 4))
	plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

	print(n)
	idxx = RV!= RV[0]
	RV = RV[idxx]
	FWHM = FWHM[idxx]
	BIS 	= BIS[idxx]
	V_span	= V_span[idxx]
	dRV_L 	= dRV_L[idxx]
	dRV_H 	= dRV_H[idxx]
	w 		= w[idxx]

	axes_1 = plt.subplot(151)
	plt.plot(RV, FWHM, 'ko', alpha=alpha)
	plt.xlabel('Jitter [m/s]')
	plt.ylabel('FWHM [km/s]')
	r1[0, n], p = stats.pearsonr(RV, FWHM)
	fit, V 	= np.polyfit(RV, FWHM, 1, w=w)
	k[0,n]  = fit			
	# delta_r = (1-r1**2)/(N-2)**0.5
	plt.title(r'$\rho = {0:.2f}$'.format(r1[0, n]), fontsize=fontsize)
	axes_1.xaxis.set_major_locator(plt.MaxNLocator(Nx))
	axes_1.yaxis.set_major_locator(plt.MaxNLocator(Ny))

	axes_2 = plt.subplot(152)
	plt.plot(RV, BIS, 'ko', alpha=alpha)
	plt.xlabel('Jitter [m/s]')
	plt.ylabel('BIS [m/s]')
	r1[1, n], p = stats.pearsonr(RV, BIS)
	fit, V 	= np.polyfit(RV, BIS, 1, w=w)
	k[1,n]  = fit		
	plt.title(r'$\rho = {0:.2f}$'.format(r1[1, n]), fontsize=fontsize)
	axes_2.xaxis.set_major_locator(plt.MaxNLocator(Nx))
	axes_2.yaxis.set_major_locator(plt.MaxNLocator(Ny))

	axes_3 = plt.subplot(153)
	plt.plot(RV, V_span, 'ko', alpha=alpha)
	plt.xlabel('Jitter [m/s]')
	plt.ylabel(r'$V_{span}$ [m/s]')
	r1[2, n], p = stats.pearsonr(RV, V_span)
	fit, V 	= np.polyfit(RV, V_span, 1, w=w)
	k[2,n]  = fit	
	plt.title(r'$\rho = {0:.2f}$'.format(r1[2, n]), fontsize=fontsize)
	axes_3.xaxis.set_major_locator(plt.MaxNLocator(Nx))
	axes_3.yaxis.set_major_locator(plt.MaxNLocator(Ny))

	axes_4 = plt.subplot(154)
	plt.plot(RV, dRV_L, 'ko', alpha=alpha)
	plt.xlabel('Jitter [m/s]')
	plt.ylabel(r'$\Delta RV_L$ [m/s]')
	r1[3, n], p = stats.pearsonr(RV, dRV_L)
	fit, V 	= np.polyfit(RV, dRV_L, 1, w=w)
	k[3,n]  = fit
	plt.title(r'$\rho = {0:.2f}; k={1:.2f}$'.format(r1[3, n],fit), fontsize=fontsize)
	axes_4.xaxis.set_major_locator(plt.MaxNLocator(Nx))
	axes_4.yaxis.set_major_locator(plt.MaxNLocator(Ny))

	axes_5 = plt.subplot(155)
	plt.plot(RV, dRV_H, 'ko', alpha=alpha)
	plt.xlabel('Jitter [m/s]')
	plt.ylabel(r'$\Delta RV_H$ [m/s]')
	r1[4, n], p = stats.pearsonr(RV, dRV_H)
	fit, V 	= np.polyfit(RV, dRV_H, 1, w=w)
	k[4,n]  = fit
	plt.title(r'$\rho = {0:.2f}; k={1:.2f}$'.format(r1[4, n],fit), fontsize=fontsize)
	axes_5.xaxis.set_major_locator(plt.MaxNLocator(Nx))
	axes_5.yaxis.set_major_locator(plt.MaxNLocator(Ny))

	plt.savefig(SAVE_NAME)
	# plt.show()
	plt.close('all')

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
		plt.plot(v_ccf, gaussian(v_ccf, *popt))
		plt.show()

		# RV, FWHM
		plt.plot(RV_HARPS, RV-RV_HARPS); plt.show()
		plt.plot(FWHM_HARPS, FWHM-FWHM_HARPS); plt.show()
		plt.plot(FWHM_HARPS, V_span); plt.show()
		
		# V_span
		plt.plot(v_CCF, gaussian(v_CCF, *popt))
		plt.plot(x_span, y_span, '.') 
		plt.show()

	#############################
	# Plot the improved version #
	#############################

if 0:

	left  	= 0.15  # the left side of the subplots of the figure
	right 	= 0.95    # the right side of the subplots of the figure
	bottom 	= 0.2   # the bottom of the subplots of the figure
	top 	= 0.8      # the top of the subplots of the figure
	wspace 	= 0.6   # the amount of width reserved for blank space between subplots
	hspace 	= 0.2   # the amount of height reserved for white space between subplots
	Nx 		= 3
	Ny 		= 3
	fontsize= 18

	plt.rcParams.update({'font.size': 20})
	fig, axes = plt.subplots(figsize=(8, 4))
	plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
	fig.suptitle('Simulation', y=0.95)

	axes_4 = plt.subplot(121)
	plt.plot(RV, dRV_L, 'ko', alpha=alpha)
	plt.xlabel('Jitter [m/s]')
	plt.ylabel(r'$\Delta RV_L$ [m/s]')
	r1, p = stats.pearsonr(RV, dRV_L)
	delta_r = (1-r1**2)/(N-2)**0.5
	fit, V 	= np.polyfit(RV, dRV_L, 1, w=w)
	plt.title(r'$\rho$={0:.2f}±{1:.2f}'.format(r1,delta_r), fontsize=fontsize)
	axes_4.xaxis.set_major_locator(plt.MaxNLocator(Nx))
	axes_4.yaxis.set_major_locator(plt.MaxNLocator(Ny))

	axes_4 = plt.subplot(122)
	plt.plot(RV, dRV_H, 'ko', alpha=alpha)
	plt.xlabel('Jitter [m/s]')
	plt.ylabel(r'$\Delta RV_H$ [m/s]')
	r1, p = stats.pearsonr(RV, dRV_H)
	delta_r = (1-r1**2)/(N-2)**0.5
	fit, V 	= np.polyfit(RV, dRV_H, 1, w=w)
	plt.title(r'$\rho$={0:.2f}±{1:.2f}'.format(r1,delta_r), fontsize=fontsize)
	axes_5.xaxis.set_major_locator(plt.MaxNLocator(Nx))
	axes_5.yaxis.set_major_locator(plt.MaxNLocator(Ny))

	plt.savefig(SAVE_NAME[:-4] + '_half_frequency.png')
	plt.show()

if 1: # size
	plt.rcParams.update({'font.size': 14})
	plt.figure()
	# plt.title('Plage')
	plt.title('Spot')
	plt.plot(S1*100, abs(r1[0, :]), 'o-', label='FWHM')
	plt.plot(S1*100, abs(r1[1, :]), 's-', label='BIS')
	plt.plot(S1*100, abs(r1[2, :]), 'p-', label=r'$V_{span}$')
	plt.plot(S1*100, abs(r1[3, :]), '^-', label=r'$\Delta RV_L$')
	plt.plot(S1*100, abs(r1[4, :]), 'v-', label=r'$\Delta RV_H$')
	plt.ylim([0,1])
	# plt.xlabel('Plage size [%]')
	plt.xlabel('Spot size [%]')
	# plt.xlabel('Spot latitude [$\degree$]')
	plt.ylabel('Pearson correlation coefficient')
	plt.grid(True)
	plt.legend(loc='right')
	plt.savefig('/Volumes/DataSSD/SOAP_2/outputs/spot_size.png')
	# plt.savefig('/Volumes/DataSSD/SOAP_2/outputs/plage_size.png')
	# plt.savefig('/Volumes/DataSSD/SOAP_2/outputs/spot_lat.png')
	plt.show()

if 1: # size
	plt.rcParams.update({'font.size': 14})
	plt.figure()
	# plt.title('Plage')
	plt.title('Spot')
	plt.plot(S1*100, abs(k[0, :]*1000), 'o-', label='FWHM')
	plt.plot(S1*100, abs(k[1, :]), 's-', label='BIS')
	plt.plot(S1*100, abs(k[2, :]), 'p-', label=r'$V_{span}$')
	plt.plot(S1*100, abs(k[3, :]), '^-', label=r'$\Delta RV_L$')
	plt.plot(S1*100, abs(k[4, :]), 'v-', label=r'$\Delta RV_H$')
	# plt.ylim([0,1])
	# plt.xlabel('Plage size [%]')
	plt.xlabel('Spot size [%]')
	plt.ylabel('k')
	plt.grid(True)
	# plt.legend(loc='best')
	plt.savefig('/Volumes/DataSSD/SOAP_2/outputs/spot_size_k.png')
	plt.show()


if 0:
	# latitude
	plt.rcParams.update({'font.size': 14})
	plt.figure()
	# plt.title('Plage')
	plt.title('Spot')
	plt.plot(lat, abs(r1[0, :]), 'o-', label='FWHM')
	plt.plot(lat, abs(r1[1, :]), 's-', label='BIS')
	plt.plot(lat, abs(r1[2, :]), 'p-', label=r'$V_{span}$')
	plt.plot(lat, abs(r1[3, :]), '^-', label=r'$\Delta RV_L$')
	plt.plot(lat, abs(r1[4, :]), 'v-', label=r'$\Delta RV_H$')
	plt.xlim([0,85])
	plt.ylim([0,1])
	plt.xlabel('Spot latitude [$\degree$]')
	# plt.xlabel('Plage latitude [$\degree$]')
	plt.ylabel('Pearson correlation coefficient')
	plt.grid(True)
	plt.legend(loc='best')
	plt.savefig('/Volumes/DataSSD/SOAP_2/outputs/spot_lat.png')
	# plt.savefig('/Volumes/DataSSD/SOAP_2/outputs/plage_lat.png')
	plt.show()

	plt.rcParams.update({'font.size': 14})
	plt.figure()
	# plt.title('Plage')
	plt.title('Spot')
	plt.plot(lat, abs(k[0, :]*1000), 'o-', label='FWHM')
	plt.plot(lat, abs(k[1, :]), 's-', label='BIS')
	plt.plot(lat, abs(k[2, :]), 'p-', label=r'$V_{span}$')
	plt.plot(lat, abs(k[3, :]), '^-', label=r'$\Delta RV_L$')
	plt.plot(lat, abs(k[4, :]), 'v-', label=r'$\Delta RV_H$')
	plt.xlim([0,85])
	plt.ylim([0,1])
	plt.xlabel('Spot latitude [$\degree$]')
	# plt.xlabel('Plage latitude [$\degree$]')
	plt.ylabel('k')
	plt.grid(True)
	plt.legend(loc='best')
	plt.savefig('/Volumes/DataSSD/SOAP_2/outputs/spot_lat_k.png')
	# plt.savefig('/Volumes/DataSSD/SOAP_2/outputs/plage_lat.png')
	plt.show()


