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
# simplify n_file @08/08/17
# Add a loop to calculate RVC and RVW @08/08/17

#############################################

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
from statistics import median

#############################################

def gaussian(x, a, mu, sigma, C):
    val = a / ( (2*math.pi)**0.5 * sigma ) * np.exp(-(x - mu)**2 / sigma**2) + C
    return val

#############################################

STAR        = 'HD216770'
# DIR         = 'HD128621_1_..2010-03-23'
# FILE        = glob.glob('../' + DIR + '/3-ccf_fits/*fits')
FILE        = glob.glob('../' + STAR + '/3-ccf_fits/*fits')
n_file      = len(FILE)
MJD         = np.zeros(n_file)
RV_g        = np.zeros(n_file)
RV_HARPS    = np.zeros(n_file)
FWHM_HARPS  = np.zeros(n_file)

for n in range(n_file):
    # progress bar #
    sys.stdout.write('\r')
    sys.stdout.write("[%-50s] %d%%" % ('='*int((n+1)*50./n_file), int((n+1)*100./n_file)))
    sys.stdout.flush()    
    
    hdulist     	= fits.open(FILE[n])
    RV_HARPS[n] 	= hdulist[0].header['HIERARCH ESO DRS CCF RVC']                 # Baryc RV (drift corrected) (km/s)
    FWHM_HARPS[n] = hdulist[0].header['HIERARCH ESO DRS CCF FWHM']

print('\n')    

RVC = median(RV_HARPS)
RVW = median(FWHM_HARPS) * 1.5
x   = np.arange(RVC-RVW, RVC+RVW+0.1, 0.1)
y   = np.zeros(len(x))
RV_noise = np.zeros(n_file)

plt.figure()

for n in range(n_file):
# for n in np.arange(n_file-3)+3:    
    
    # progress bar #
    sys.stdout.write('\r')
    sys.stdout.write("[%-50s] %d%%" % ('='*int((n+1)*50./n_file), int((n+1)*100./n_file)))
    sys.stdout.flush()    
    
    hdulist     = fits.open(FILE[n])
    v0          = hdulist[0].header['CRVAL1']                                   # velocity on the left (N_STARting point)
    print(v0)
    STAR_read   = hdulist[0].header['OBJECT']
    
    # remove file if necessary
    # if (STAR_read != STAR) and (STAR_read != 'HD-10700') and (STAR_read != 'HR509') and (STAR_read != 'TauCeti') \
    #     and (STAR_read != 'Tau-Cet') and (STAR_read != 'TAUCET') and (STAR_read != 'tau-Cet') \
    #     and (STAR_read != 'Tau-Ceti') and (STAR_read != 'HD049933'):
    #     print (' Achtung! ' + STAR_read + ' instead of '+ STAR)
    #     shutil.move(FILE[n], '../' + STAR + '/3-ccf_fits/abandoned/')
    #     continue


    # if (STAR_read != STAR) and (STAR_read != 'HR5460'):
    #     print (' Achtung! ' + STAR_read + ' instead of '+ STAR)
    #     shutil.move(FILE[n], '../' + STAR + '/3-ccf_fits/abandoned/')
    #     continue        

    RV_HARPS[n] = hdulist[0].header['HIERARCH ESO DRS CCF RVC']                 # Baryc RV (drift corrected) (km/s)
    if (RV_HARPS[n] > RVC+10) or (RV_HARPS[n] < RVC-10):
        print(' Achtung! ' +  STAR_read + ' largely offset')
        shutil.move(FILE[n], '../' + STAR + '/3-ccf_fits/abandoned/')
        continue

    RV_noise[n] = hdulist[0].header['HIERARCH ESO DRS CCF NOISE'] * 1000        # RV_noise in m/s    
    if RV_noise[n] > 5:
        print(' Achtung! ' + STAR_read + ' too noisy')
        shutil.move(FILE[n], '../' + STAR + '/3-ccf_fits/abandoned/')
        continue

    # if (v0 > 11.1) or (v0 < 10.9):
    #     print(' Achtung! ' +  STAR_read + ' initial velocity mismatch ')
    #     shutil.move(FILE[n], '../' + STAR + '/3-ccf_fits/abandoned/')
    #     continue

    if 1: # test
        quartile_1, quartile_3 = np.percentile(RV_HARPS, [25, 75])
        iqr = quartile_3 - quartile_1
        lower_bound = quartile_1 - (iqr * 10)
        upper_bound = quartile_3 + (iqr * 10)
        if (RV_HARPS[n] < lower_bound) or (RV_HARPS[n] > upper_bound):
            print(' Achtung! ' +  STAR_read + ' odd offset')
            shutil.move(FILE[n], '../' + STAR + '/3-ccf_fits/abandoned/')
            continue

    CCF         = hdulist[0].data                                               # ccf 2-d array
    ccf         = CCF[- 1, :]                                                   # ccf 1-d array (whole range)
    delta_v     = hdulist[0].header['CDELT1']                                   # velocity grid size 
    v           = v0 + np.arange(CCF.shape[1]) * delta_v                        # velocity array (whole range)
    # plt.plot(v,ccf)
    print(min(ccf)**0.5)

    f           = CubicSpline(v, ccf / ccf.max())
    y           = f(x)
    popt, pcov  = curve_fit( gaussian, x, y, [-RVW/2, RV_HARPS[n], RVW/2, 1])
    # plt.plot(x, y, x, gaussian(x, *popt))
    RV_g[n]     = popt[1]
    
    if abs(RV_HARPS[n] - RV_g[n])*1000 > 5:
        print(' Achtung! ' +  STAR_read + ' RV fitting issue')
        shutil.move(FILE[n], '../' + STAR + '/3-ccf_fits/abandoned/')    
        continue

    if 0:
    # valid for HD128621_3_2010-06-13..2012-12-31
        if RV_HARPS[n] *1000 > 250:
    # valid for HD128621_1_..2010-03-23
        # if (RV_HARPS[n] - np.mean(RV_HARPS)) *1000 > 100:
            print(' Achtung! ' +  STAR_read + ' outlier')
            shutil.move(FILE[n], '../' + STAR + '/3-ccf_fits/abandoned/')    
            continue            

    x_new       = x
    y_new       = (y - popt[3]) / popt[0]

    if 0:
        # set valid to HD22049 (epslon Eri)
        # if max(y_new>0.12): 
        # valid for HD128621_1_..2010-03-23
        if max(y_new>0.11):
            print(' Achtung! ' +  STAR_read + ' outlier?')
            shutil.move(FILE[n], '../' + STAR + '/3-ccf_fits/abandoned/')    
            continue    

    MJD[n]      = hdulist[0].header['MJD-OBS']                              # modified Julian Date (JD - 2400000.5)

    if 0: 
    # valid for HD189733
        if (MJD[n] < 53986.01) or (MJD[n] > 53990):
    # valid for HD128621_3_2010-06-13..2012-12-31
        # if ((MJD[n] > 55700) or (MJD[n] < 55600)):
    # valid for HD128621_1_..2010-03-23
        # if ((MJD[n] > 55100) or (MJD[n] < 54850)):        
            print(' Achtung! ' +  STAR_read + ' not selected')
            shutil.move(FILE[n], '../' + STAR + '/3-ccf_fits/abandoned/')    
            continue

    if 0:
        if (MJD[n] > 57161):
            print(' Achtung! ' +  STAR_read + ' fibre upgrade')
            shutil.move(FILE[n], '../' + STAR + '/3-ccf_fits/abandoned/')    
            continue

    plt.plot(x_new, y_new, '-')
    output_name = FILE[n]
    output_name = output_name.replace('3-ccf_fits', '4-ccf_dat')
    writefile   = output_name.replace('.fits', '.dat')
    np.savetxt(writefile, y_new)

plt.show()
print('\n')        
idx = (MJD == 0)
# idx = (MJD > 55100) | (MJD < 54850) | (MJD == 0);
np.savetxt('../' + STAR + '/info.dat', [RVC, RVW])
np.savetxt('../' + STAR + '/MJD.dat', MJD[~idx])
np.savetxt('../' + STAR + '/RV_HARPS.dat', RV_HARPS[~idx])
np.savetxt('../' + STAR + '/x.dat', x)
np.savetxt('../' + STAR + '/RV_noise.dat', RV_noise[~idx])
np.savetxt('../' + STAR + '/FWHM.dat', FWHM_HARPS[~idx])
output  = np.vstack((MJD[~idx], (RV_HARPS[~idx]-np.mean(RV_HARPS[~idx]))*1000, RV_noise[~idx]))
output = np.transpose(output)
np.savetxt('../' + STAR + '/' + STAR + '.txt', output, fmt='%1.8f')


# Verification # 
if 0:
    plt.errorbar(MJD[~idx], RV_g[~idx] *1000 - np.mean(RV_g[~idx] *1000), yerr=RV_noise[~idx], fmt=".k", capsize=0)
    plt.errorbar(MJD, RV_HARPS *1000, yerr=RV_noise, fmt=".k", capsize=0)
    plt.errorbar(MJD[~idx], RV_HARPS[~idx] *1000 - np.mean(RV_HARPS[~idx] *1000), yerr=RV_noise[~idx], fmt=".k", capsize=0, alpha=0.3)
    plt.show()
    # idx2 = (MJD[~idx] > 53986) & (MJD[~idx] < 53990)
    # plt.errorbar(MJD[~idx][idx2], RV_HARPS[~idx][idx2] *1000 - np.mean(RV_HARPS[~idx][idx2] *1000), yerr=RV_noise[~idx][idx2], fmt=".k", capsize=0, alpha=0.3)
    # plt.show()

if 0:
    from numpy.polynomial import polynomial as P
    c, stats    = P.polyfit(MJD[~idx], RV_HARPS[~idx] *1000 - np.mean(RV_HARPS[~idx] *1000) ,3, full=True, w = 1/(RV_noise[~idx])**2)
    x_fit       = np.linspace(min(MJD[~idx]-1), max(MJD[~idx]+1), 10000)
    y_fit       = P.polyval(x_fit, c)
    plt.errorbar(MJD[~idx], RV_HARPS[~idx] *1000 - np.mean(RV_HARPS[~idx] *1000), yerr=RV_noise[~idx], fmt=".k", capsize=0)
    plt.plot(x_fit, y_fit)
    plt.show()

    y_plot = P.polyval(MJD[~idx], c)
    plt.plot(MJD[~idx], y_plot - (RV_HARPS[~idx] *1000 - np.mean(RV_HARPS[~idx] *1000)), '.')
    plt.show()








