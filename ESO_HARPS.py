#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 16:43:04 2017

@author: jzhao
"""

"""
    Extract CCF: .fits -> .dat
"""

###########
# Updates #
###########
# Remove hdulist[0].header['HIERARCH ESO DRS BERV']. 

# import os
import sys
import shutil
from astropy.io import fits
import glob
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
import math

def gaussian(x, a, mu, sigma, C):
    val = a / ( (2*math.pi)**0.5 * sigma ) * np.exp(-(x - mu)**2 / sigma**2) + C
    return val

# Root mean square
#def rms(num):
#    return sqrt(sum(n*n for n in num)/len(num))
#############################################

star        = 'Gl581'
FILE        = glob.glob('../' + star + '/3-ccf_fits/*fits')
N           = len(FILE)
N_start     = 0
N_end       = N
n_file      = N_end - N_start
MJD         = np.zeros(n_file)
RV_g        = np.zeros(n_file)
RV_HARPS    = np.zeros(n_file)

# x           = np.arange(12.5-10, 12.5+10+0.1, 0.1)                            # HD96700 FWHM = 6.9      # over sampling to 0.1 km/s [-10.2, -0.8]
# x           = np.arange(5.5-18, 5.5+18+0.1, 0.1)                                # HD142 FWHM = 15
# x           = np.arange(31-10, 31+10+0.1, 0.1)                                  # CoRot-7
x           = np.arange(-9.2-4.5, -9.2+4.5+0.1, 0.1) 
y           = np.zeros(len(x))

plt.figure()

for n in range(N_start, N_end):
    
    # progress bar #
    sys.stdout.write('\r')
    sys.stdout.write("[%-50s] %d%%" % ('='*int((n+1-N_start)*50./(N_end-N_start)), int((n+1-N_start)*100./(N_end-N_start))))
    sys.stdout.flush()    
    
    hdulist     = fits.open(FILE[n])
    v0          = hdulist[0].header['CRVAL1']                                   # velocity on the left (N_starting point)
    RV_HARPS[n] = hdulist[0].header['HIERARCH ESO DRS CCF RVC']                 # Baryc RV (drift corrected) (km/s)
    v_noise     = hdulist[0].header['HIERARCH ESO DRS CCF NOISE'] * 1000        # RV_noise in m/s
    star_read   = hdulist[0].header['OBJECT']
    
    # remove file if necessary
    if star_read != star:
        print (' Achtung! ' + star_read + ' instead of '+ star)
        shutil.move(FILE[n], '../' + star + '/3-ccf_fits/abandoned/')
        continue
    
    if v_noise > 5:
        print(' Achtung! ' + star_read + ' too noisy')
        shutil.move(FILE[n], '../' + star + '/3-ccf_fits/abandoned/')
        continue
    
    if (RV_HARPS[n] > -9.2+10) or (RV_HARPS[n] < -9.2-10):
        print(' Achtung! ' +  star_read + ' largely offset')
        shutil.move(FILE[n], '../' + star + '/3-ccf_fits/abandoned/')
        continue

    MJD[n]      = hdulist[0].header['MJD-OBS']
    
    CCF         = hdulist[0].data                                               # ccf 2-d array
    ccf         = CCF[- 1, :]                                                   # ccf 1-d array (whole range)
    delta_v     = hdulist[0].header['CDELT1']                                   # velocity grid size 
    v           = v0 + np.arange(CCF.shape[1]) * delta_v                        # velocity array (whole range)
    # plt.plot(v,ccf)

    f           = CubicSpline(v, ccf / ccf.max())
    y           = f(x)
    # popt, pcov  = curve_fit( gaussian, x, y, [-1.76, RV_HARPS[n], 2.5, 1])
    popt, pcov  = curve_fit( gaussian, x, y, [-5, RV_HARPS[n], 4, 1])
    # plt.plot(x, y, x, gaussian(x, *popt))
    RV_g[n]     = popt[1]
    
    if abs(RV_HARPS[n] - RV_g[n])*1000 > 5:
        print(' Achtung! ' +  star_read + ' RV fitting issue')
        shutil.move(FILE[n], '../' + star + '/3-ccf_fits/abandoned/')    
        continue
        
    x_new       = x
    y_new       = (y - popt[3]) / popt[0]
    plt.plot(x_new, y_new, '-')
    
    output_name = FILE[n]
    output_name = output_name.replace('3-ccf_fits', '4-ccf_dat')
    writefile   = output_name.replace('.fits', '.dat')
    np.savetxt(writefile, y_new)

# Verification # 
if 0:
    idx = (RV_g == 0)
    plt.plot( RV_g[~idx] *1000, RV_HARPS[~idx] * 1000, '+')
    plt.plot(MJD[~idx], RV_g[~idx] *1000, '.', MJD[~idx], RV_HARPS[~idx] * 1000, '+')
    plt.plot(MJD, (RV_g - np.mean(RV_g))*1000, '.', MJD, (RV_HARPS - np.mean(RV_HARPS)) * 1000, '+')
    plt.show()