#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 10:46:40 2017

@author: jzhao
"""



import sys
import shutil
from astropy.io import fits
import glob
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
import math
from statistics import median

def gaussian(x, a, mu, sigma, C):
    val = a / ( (2*math.pi)**0.5 * sigma ) * np.exp(-(x - mu)**2 / sigma**2) + C
    return val

# Root mean square
#def rms(num):
#    return sqrt(sum(n*n for n in num)/len(num))
#############################################


STAR        = 'HD189733'
FILE        = glob.glob('../' + STAR + '/3-ccf_fits/*fits')
n_file      = len(FILE)
MJD         = np.zeros(n_file)
RV_g        = np.zeros(n_file)
RV_HARPS    = np.zeros(n_file)
FWHM_HARPS  = np.zeros(n_file)
ccf_max     = np.zeros(n_file)

for n in range(n_file):
    # progress bar #
    sys.stdout.write('\r')
    sys.stdout.write("[%-50s] %d%%" % ('='*int((n+1)*50./n_file), int((n+1)*100./n_file)))
    sys.stdout.flush()    
    
    hdulist     	= fits.open(FILE[n])
    RV_HARPS[n] 	= hdulist[0].header['HIERARCH ESO DRS CCF RVC']                 # Baryc RV (drift corrected) (km/s)
    FWHM_HARPS[n] = hdulist[0].header['HIERARCH ESO DRS CCF FWHM']

print('\n')  

RVC     = median(RV_HARPS)
RVW     = median(FWHM_HARPS) * 1.5
x       = np.arange(RVC-RVW, RVC+RVW+0.1, 0.1)
y       = np.zeros(len(x))
x_tmp   = np.arange(RVC-RVW, RVC+RVW+0.1, 0.1)
Y_tmp   = np.arange(RVC-RVW, RVC+RVW+0.1, 0.1)

plt.figure()

for n in range(n_file):
    
    # progress bar #
    sys.stdout.write('\r')
    sys.stdout.write("[%-50s] %d%%" % ('='*int((n+1)*50./n_file), int((n+1)*100./n_file)))
    sys.stdout.flush()    
    
    hdulist     = fits.open(FILE[n])
    v0          = hdulist[0].header['CRVAL1']                                   # velocity on the left (N_STARting point)
    RV_HARPS[n] = hdulist[0].header['HIERARCH ESO DRS CCF RVC']                 # Baryc RV (drift corrected) (km/s)
    v_noise     = hdulist[0].header['HIERARCH ESO DRS CCF NOISE'] * 1000        # RV_noise in m/s
    STAR_read   = hdulist[0].header['OBJECT']
    
    # remove file if necessary
    if STAR_read != STAR:
        print (' Achtung! ' + STAR_read + ' instead of '+ STAR)
        shutil.move(FILE[n], '../' + STAR + '/3-ccf_fits/abandoned/')
        continue
    
    if v_noise > 5:
        print(' Achtung! ' + STAR_read + ' too noisy')
        shutil.move(FILE[n], '../' + STAR + '/3-ccf_fits/abandoned/')
        continue
    
    if (RV_HARPS[n] > -9.2+10) or (RV_HARPS[n] < -9.2-10):
        print(' Achtung! ' +  STAR_read + ' largely offset')
        shutil.move(FILE[n], '../' + STAR + '/3-ccf_fits/abandoned/')
        continue

    MJD[n]      = hdulist[0].header['MJD-OBS']
    
    CCF         = hdulist[0].data                                               # ccf 2-d array
    ccf         = CCF[- 1, :]                                                   # ccf 1-d array (whole range)
    delta_v     = hdulist[0].header['CDELT1']                                   # velocity grid size 
    v           = v0 + np.arange(CCF.shape[1]) * delta_v                        # velocity array (whole range)
    # plt.plot(v,ccf)

    f           = CubicSpline(v, ccf / ccf.max())
    y           = f(x)
    popt, pcov  = curve_fit( gaussian, x, y, [-1.8, RV_HARPS[n], 2.5, 1])
    # plt.plot(x, y, x, gaussian(x, *popt))
    RV_g[n]     = popt[1]
    
    if abs(RV_HARPS[n] - RV_g[n])*1000 > 5:
        print(' Achtung! ' +  STAR_read + ' RV fitting issue')
        shutil.move(FILE[n], '../' + STAR + '/3-ccf_fits/abandoned/')    
        continue
        
    x_new       = x - popt[1]
    y_new       = (y - popt[3]) / popt[0]
    plt.plot(x_new, y_new, '-')
    ccf_max[n]  = ccf.max()
    
    f           = CubicSpline(x_new, y_new)
    y_tmp       = f(x_tmp) *  ccf.max()
    Y_tmp       = Y_tmp + y_tmp
#    plt.plot(x_tmp, f(x_tmp) - Y_tmp)
    
Y_tmp = Y_tmp / sum(ccf_max)

popt, pcov  = curve_fit( gaussian, x_tmp, Y_tmp)

writefile = ('../' + STAR + '/template1.dat')
np.savetxt(writefile, Y_tmp)

#    output_name = FILE[n]
#    output_name = output_name.replace('3-ccf_fits', '4-ccf_dat')
#    writefile   = output_name.replace('.fits', '.dat')
#    np.savetxt(writefile, y_new)

# Verification # 
if 0:
    idx = (RV_g == 0)
    plt.plot( RV_g[~idx] *1000, RV_HARPS[~idx] * 1000, '+')
    plt.plot(MJD[~idx], RV_g[~idx] *1000, '.', MJD[~idx], RV_HARPS[~idx] * 1000, '+')
    plt.plot(MJD, (RV_g - np.mean(RV_g))*1000, '.', MJD, (RV_HARPS - np.mean(RV_HARPS)) * 1000, '+')
    plt.show()