#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 10:46:40 2017

@author: jzhao
"""
# Update [-RVW/2, RV_HARPS[n], RVW/2, 1] 22/08/17


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

#############################################
def gaussian(x, a, mu, sigma, C):
    val = a / ( (2*math.pi)**0.5 * sigma ) * np.exp(-(x - mu)**2 / sigma**2) + C
    return val
#############################################

STAR        = 'HD189733'
FILE        = glob.glob('../' + STAR + '/3-ccf_fits/*fits')
n_file      = len(FILE)
MJD         = np.zeros(n_file)
RV_g        = np.zeros(n_file)
RV_HARPS    = np.zeros(n_file)
FWHM_HARPS  = np.zeros(n_file)
y_max       = np.zeros(n_file)

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
x       = np.arange(RVC-RVW, RVC+RVW, 0.1)
y       = np.zeros(len(x))
x_tmp   = np.arange(RVC-RVW, RVC+RVW, 0.1)
Y_tmp   = np.zeros(len(x_tmp))

plt.figure()

for n in range(n_file):
    
    # progress bar #
    sys.stdout.write('\r')
    sys.stdout.write("[%-50s] %d%%" % ('='*int((n+1)*50./n_file), int((n+1)*100./n_file)))
    sys.stdout.flush()    
    
    hdulist     = fits.open(FILE[n])
    v0          = hdulist[0].header['CRVAL1']                                   # velocity on the left (N_STARting point)
    v_noise     = hdulist[0].header['HIERARCH ESO DRS CCF NOISE'] * 1000        # RV_noise in m/s
    STAR_read   = hdulist[0].header['OBJECT']
    MJD[n]      = hdulist[0].header['MJD-OBS']
    
    # remove file if necessary
    if STAR_read != STAR:
        print (' Achtung! ' + STAR_read + ' instead of '+ STAR)
        shutil.move(FILE[n], '../' + STAR + '/3-ccf_fits/abandoned/')
        continue
    
    if v_noise > 5:
        print(' Achtung! ' + STAR_read + ' too noisy')
        shutil.move(FILE[n], '../' + STAR + '/3-ccf_fits/abandoned/')
        continue
    
    if (RV_HARPS[n] > RVC+10) or (RV_HARPS[n] < RVC-10):
        print(' Achtung! ' +  STAR_read + ' largely offset')
        shutil.move(FILE[n], '../' + STAR + '/3-ccf_fits/abandoned/')
        continue

    CCF         = hdulist[0].data                                               # ccf 2-d array
    ccf         = CCF[- 1, :]                                                   # ccf 1-d array (whole range)
    delta_v     = hdulist[0].header['CDELT1']                                   # velocity grid size 
    
    v           = v0 + np.arange(CCF.shape[1]) * delta_v                        # velocity array (whole range)
    idx_v       = (v > RVC-RVW-0.1) & (v < RVC+RVW+0.1)
    x           = v[idx_v]
    y           = ccf[idx_v]
    y_max[n]    = y.max()
    
    popt, pcov  = curve_fit( gaussian, x, y/max(y), [-RVW/2, RV_HARPS[n], RVW/2, 1])
    RV_g[n]     = popt[1]
    
    if abs(RV_HARPS[n] - RV_g[n])*1000 > 5:
        print(' Achtung! ' +  STAR_read + ' RV fitting issue')
        shutil.move(FILE[n], '../' + STAR + '/3-ccf_fits/abandoned/')    
        continue

    f           = CubicSpline( x-(popt[1]-RVC), y )                             # shift the observed spectrum in order to co-add to a template
    y_tmp       = f(x_tmp)
    Y_tmp       = Y_tmp + y_tmp
#    plt.figure()
#    plt.plot(x_tmp - RVC, y_tmp / y.max() - Y_tmp)
#    plt.plot(x_tmp - RVC, y_tmp)
    
Y_tmp = Y_tmp / sum(y_max)

writefile = ('../' + STAR + '/template1.dat')
np.savetxt(writefile, Y_tmp)


# Verification # 
if 0:
    idx = (RV_g == 0)
    plt.plot( RV_g[~idx] *1000, RV_HARPS[~idx] * 1000, '+')
    plt.plot(MJD[~idx], RV_g[~idx] *1000, '.', MJD[~idx], RV_HARPS[~idx] * 1000, '+')
    plt.plot(MJD, (RV_g - np.mean(RV_g))*1000, '.', MJD, (RV_HARPS - np.mean(RV_HARPS)) * 1000, '+')
    plt.show()