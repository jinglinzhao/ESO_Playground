#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 15:01:51 2017

@author: jzhao
"""
# modify to find out the 1st order coefficient @25/08/17
# further modify to find out the shift of the observation so that the 1st order coefficient  becomes 0 @30/08/17
# introduce error propropagation and weights @30/08/17


import sys
#import shutil
from astropy.io import fits
import glob
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
import math
#from math import  exp
from statistics import median
from scipy.optimize import minimize
from numpy.polynomial.hermite import hermfit

#############################################
def gaussian(x, a, mu, sigma, C):
    val = a / ( (2*math.pi)**0.5 * sigma ) * np.exp(-(x - mu)**2 / sigma**2) + C
    return val

# Root mean square
def rms(num):
    return math.sqrt(sum(n*n for n in num)/len(num))

#############################################

STAR        = 'HD189733'
FILE        = glob.glob('../' + STAR + '/3-ccf_fits/*fits')
n_file      = len(FILE)
MJD, RV_fit, RV_HARPS, FWHM_HARPS, y_max    = np.zeros((5, n_file))
LPD, A, S, K                                = np.zeros((4, n_file))
array1      = np.zeros(n_file)

print('Reading data...')
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
RVW     = median(FWHM_HARPS) * 1.4
x       = np.arange(RVC-RVW, RVC+RVW, 0.1)
y       = np.zeros(len(x))
x_tmp   = np.arange(RVC-RVW, RVC+RVW, 0.1)
Y_tmp1  = np.loadtxt('../' + STAR + '/template1.dat')
Y_tmp1_err = np.loadtxt('../' + STAR + '/template1_err.dat')
Y_tmp1_nor = Y_tmp1 / max(Y_tmp1)
Y_tmp1_err = Y_tmp1_err / max(Y_tmp1)
#Y_tmp2  = np.zeros(len(x_tmp))

#plt.figure()
#Y_tmp1_err  = Y_tmp1**0.5
popt, pcov  = curve_fit( gaussian, x_tmp, Y_tmp1_nor, [-RVW/2, RVC, RVW/2, 1])
Y_tmp1_nor  = (Y_tmp1_nor - popt[3]) / popt[0]
Y_tmp1_err  = Y_tmp1_err / abs(popt[0])

print('Calculating...')
for n in range(n_file):
    # progress bar #
    sys.stdout.write('\r')
    sys.stdout.write("[%-50s] %d%%" % ('='*int((n+1)*50./n_file), int((n+1)*100./n_file)))
    sys.stdout.flush()    
    
    hdulist     = fits.open(FILE[n])
    v0          = hdulist[0].header['CRVAL1']                                   # velocity on the left (N_STARting point)
    MJD[n]      = hdulist[0].header['MJD-OBS']
    
    CCF         = hdulist[0].data                                               # ccf 2-d array
    ccf         = CCF[- 1, :]                                                   # ccf 1-d array (whole range)
    delta_v     = hdulist[0].header['CDELT1']                                   # velocity grid size 
    
    v           = v0 + np.arange(CCF.shape[1]) * delta_v                        # velocity array (whole range)
    idx_v       = (v > RVC-RVW-0.1) & (v < RVC+RVW+0.1)
    x           = v[idx_v]
    y           = ccf[idx_v] 
    y_max[n]    = y.max()

    # test begins
    def obj(x_shift):
        f_test      = CubicSpline(v + x_shift, ccf)
        y_test      = f_test(x_tmp)
        y_test_err  = y_test**0.5
        y_test      = y_test / max(y_test)
        y_test_err  = y_test_err / max(y_test)
        y_test      = (y_test - popt[3]) / popt[0]
        y_test_err  = y_test_err / abs(popt[0])
    #    plt.figure(); plt.plot(x_tmp, y_test, x_tmp, Y_tmp1)
        array_hermite = hermfit(x_tmp - popt[1], y_test / Y_tmp1_nor, 4, w = 1 / (y_test/Y_tmp1_nor * ((y_test/y_test_err)**2 + (Y_tmp1_err/Y_tmp1_nor)**2)**0.5)**2)
        return abs(array_hermite[1])
    res = minimize(obj, 0)
    array1[n] = res.x[0]
    
plt.figure();
plt.plot(MJD, array1, '.')
plt.xlabel('MJD')
plt.ylabel('x_shift4')



    
#    y           = ccf[idx_v] / max(ccf[idx_v])
#    popt, pcov  = curve_fit( gaussian, x, y, [-RVW/2, RV_HARPS[n], RVW/2, 1])
#    y           = (y - popt[3]) / popt[0]                                       # x is in original v grid; y is normalized

        
    # Now shift the template to match the observation
'''
    # Option 1: minimize rms
    def obj(x_shift):
        f           = CubicSpline(x_tmp + x_shift, Y_tmp1)                      # shift the template to the right (for x_shift > 0) to match the observation
        y_shift     = f(x)
        return rms( (y_shift * max(ccf[idx_v]) - y) / ccf[idx_v] )
    res     = minimize(obj, 0) 
    RV_fit[n] = res.x[0]
    
    f           = CubicSpline( (x - res.x[0]), y )                              # shift the observed spectrum to the left (for x_shift > 0) in order to co-add to the template
    y_tmp       = f(x_tmp)
#    Y_tmp2       = Y_tmp2 + y_tmp
    
#    plt.figure()
#    plt.plot(x_tmp - RVC, y_tmp / y.max() - Y_tmp2)
    
    obs         = y_tmp / y.max()
    x_new       = x_tmp - RVC
    LPD[n]      = sum((y_tmp / y.max() - Y_tmp2) * x_new)
    sigma_2     = sum(obs * x_new**2) / sum(obs)
    A[n]        = 1. / sigma_2**0.5 * sum(obs * x_new**1) / sum(obs)
    S[n]        = 1. / sigma_2**1.5 * sum(obs * x_new**3) / sum(obs)
    K[n]        = 1. / sigma_2**2.0 * sum(obs * x_new**4) / sum(obs)   
'''
    # Option 2: maximize exponential
#    def obj(x_shift):
#        f           = CubicSpline(x_tmp + x_shift, Y_tmp1)
#        y_shift     = f(x)
#        return -sum(exp( -((y_shift - y)/0.001)**2 ))
#    
#    res     = minimize(obj, 1)
#    RV_fit[n] = res.x[0]
        
#plt.plot(MJD, RV_fit, '.', MJD, RV_HARPS-RVC , '+')   
#Y_tmp2 = Y_tmp2 / sum(y_max)
    



'''    
    
    RV_g[n]     = popt[1]

    x_new       = x - (popt[1] - RVC)
    y_new       = (y - popt[3]) / popt[0]
    ccf_max[n]  = ccf.max()
    
    f           = CubicSpline(x_new, y_new)
    y_tmp       = f(x_tmp)
    corr = np.correlate(y_tmp, Y_tmp1, 'same')
    popt, pcov  = curve_fit( gaussian, x_tmp, corr/max(corr))
    
    x_new       = x_new + popt[1]
    print(popt[1]*1000)
    f           = CubicSpline(x_new, y_new)
    y_tmp       = f(x_tmp) * ccf_max[n]
    
    Y_tmp2       = Y_tmp2 + y_tmp
    
Y_tmp2 = Y_tmp2 / sum(ccf_max)
plt.plot(x_tmp, Y_tmp2)
popt, pcov  = curve_fit( gaussian, x_tmp, Y_tmp2)

writefile = ('../' + STAR + '/template2.dat')
np.savetxt(writefile, Y_tmp2)
'''