#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 15:01:51 2017

@author: jzhao
"""
# modify to find out the 1st order coefficient @25/08/17
# further modify to find out the shift of the observation so that the 1st order coefficient  becomes 0 @30/08/17
# introduce error propropagation and weights @30/08/17
# correct normalization error -> y_test_err  = y_test**0.5 / max(y_test) @31/08/17

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
from numpy.polynomial.hermite import hermval
from scipy.optimize import curve_fit

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
array2      = np.zeros(n_file)
array3      = np.zeros(n_file)
array_popt  = np.zeros((6, n_file))

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
RVW     = median(FWHM_HARPS) * 1.3
buffer  = max(RV_HARPS) - min(RV_HARPS)
x_tmp   = np.arange(-RVW, RVW+0.1, 0.1)
Y_tmp1  = np.loadtxt('../' + STAR + '/template1.dat')
Y_tmp1_nor = Y_tmp1 / max(Y_tmp1)
Y_tmp1_err = Y_tmp1**0.5 / max(Y_tmp1)
#Y_tmp2  = np.zeros(len(x_tmp))

#plt.figure()
popt, pcov  = curve_fit( gaussian, x_tmp, Y_tmp1_nor, [-RVW/2, RVC, RVW/2, 1], sigma = Y_tmp1_err)
Y_tmp1_nor  = (Y_tmp1_nor - popt[3]) / popt[0] * popt[2]
Y_tmp1_err  = Y_tmp1_err / abs(popt[0]) * abs(popt[2])

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
    idx_v       = (v > RVC - RVW - buffer) & (v < RVC + RVW + delta_v + buffer)
    x           = v[idx_v]
    y           = ccf[idx_v] 
    y_max[n]    = y.max()

    '''
    # test1: find x_shift to make 1st order coefficient zero
    def obj(x_shift):
        f_test      = CubicSpline(v - RVC - x_shift, ccf)
        y_test      = f_test(x_tmp)
        y_test_err  = y_test**0.5 / max(y_test)                                 # Achtung!
        y_test      = y_test / max(y_test)
        y_test      = (y_test - popt[3]) / popt[0]                              # same normalization with the template
        y_test_err  = y_test_err / abs(popt[0])
    #    plt.figure(); plt.plot(x_tmp, y_test, x_tmp, Y_tmp1)
        array_hermite = hermfit(x_tmp / popt[2], y_test / Y_tmp1_nor, 5, \
                                w = 1 / (y_test/Y_tmp1_nor * ((y_test_err/y_test)**2 + (Y_tmp1_err/Y_tmp1_nor)**2)**0.5)**2)
        return abs(array_hermite[1])
    res = minimize(obj, 0)
    array1[n] = res.x[0]
    '''
    
    # test3: curve fit
    f_obs   = CubicSpline(v - RVC, ccf)
    y_obs   = f_obs(x_tmp)
    y_obs_err   = y_obs**0.5 / max(y_obs) 
    y_obs   = y_obs / max(y_obs)
    y_obs   = (y_obs - popt[3]) / popt[0] * popt[2]
    y_obs_err   = y_obs_err / abs(popt[0]) * abs(popt[2])
    
    def obj3(x, x_shift, a0, a2, a3, a4, a5):
        model   = gaussian( (x-x_shift)/popt[2], 1, 0, 1, 0 ) * hermval((x-x_shift)/popt[2], [a0, 0, a2, a3, a4, a5])
        return model
#    model   = gaussian( (x_tmp-x_shift)/popt[2], 1, 0, 1, 0 ) * hermval((x_tmp-x_shift)/popt[2], [a0, 0, a2, a3, a4, a5])
#    w = 1 / (y_obs/Y_tmp1_nor * ((y_obs_err/y_obs)**2 + (Y_tmp1_err/Y_tmp1_nor)**2)**0.5)**2
#    resd    = (y_obs - model)**2 * w
    popt3, pcov3 = curve_fit(obj3, x_tmp, y_obs, sigma=y_obs_err, bounds=(-0.5, [0.5, 1.5, 0.5, 0.5, 0.5, 0.5]))
#    plt.figure(); plt.plot(x_tmp, y_obs, 'b', x_tmp, obj3(x_tmp, *popt3), 'r--')
    array_popt[:,n] = popt3
#    res3 = minimize(obj3, 0, 0, 0, 0, 0, 0, method='BFGS')


#==============================================================================
# Output figures
#==============================================================================

plt.figure()
plt.plot(MJD, (array_popt[0,:] )*1000, '.', MJD , (RV_HARPS-RVC)*1000, '^')
#plt.plot(range(N), (array_popt[0,:]-v_random)*1000, '.', range(N), (rv_g-v_random)*1000, '^')
plt.xlabel('MJD')
plt.ylabel('RV [m/s]')

plt.figure()
plt.plot( (RV_HARPS- np.mean(RV_HARPS))*1000, (array_popt[0,:] - np.mean(array_popt[0,:]))*1000, '.')
# plt.plot((rv_g-v_random)*1000, (array_popt[0,:]-v_random)*1000, '.')
plt.xlabel('HARPS RV [m/s]')
plt.ylabel('G-H RV [m/s]')

for i in range(1,6):
    plt.figure()
    plt.plot(MJD, array_popt[i,:], '.')
    plt.xlabel('MJD')
    plt.ylabel('popt')
    if i == 1:
        plt.title('0')    
    else:
        plt.title(str(i))
    
    
#plt.figure();
#plt.plot(MJD, RV_HARPS - np.mean(RV_HARPS), '.', MJD, array3 - np.mean(array3), 'o', MJD, array1 - np.mean(array1), '^')
#plt.xlabel('MJD')
#plt.ylabel('RV')
#plt.title('Fitting up to 5th order')    

    
'''    @ 22/09/17
    
    if 0:
        # test2: curve fit; not working 
        def obj2(x_shift, a0, a2, a3, a4, a5):
            f_test      = CubicSpline(v - RVC - x_shift, ccf)
            y_test      = f_test(x_tmp)
            y_test_err  = y_test**0.5 / max(y_test)                                 # Achtung!
            y_test      = y_test / max(y_test)
            y_test      = (y_test - popt[3]) / popt[0]                              # same normalization with the template
            y_test_err  = y_test_err / abs(popt[0])
            residual    = (y_test - Y_tmp1_nor * hermval(x_tmp / popt[2], [a0, 0, a2, a3, a4, a5])) \
                            / (y_test/Y_tmp1_nor * ((y_test_err/y_test)**2 + (Y_tmp1_err/Y_tmp1_nor)**2)**0.5)**2
            return residual
        res2 = minimize(obj2, 0, 0, 0, 0, 0, 0, method='BFGS', tol=1e-6)
        array2[n] = res2.x[0]
    #        popt2, pcov2 = curve_fit(func, xdata, ydata, bounds=(0, [3., 2., 1.]))
    
    
    
    
    
    
    
    f_test      = CubicSpline(v - RVC - res.x[0], ccf)
    y_test      = f_test(x_tmp)
    y_test_err  = y_test**0.5 / max(y_test)                                 # Achtung!
    y_test      = y_test / max(y_test)
    y_test      = (y_test - popt[3]) / popt[0]                              # same normalization with the template
    y_test_err  = y_test_err / abs(popt[0])    
    array_hermite = hermfit(x_tmp, y_test / Y_tmp1_nor, 5, w = 1 / (y_test/Y_tmp1_nor * ((y_test_err/y_test)**2 + (Y_tmp1_err/Y_tmp1_nor)**2)**0.5)**2)
    
    y_corrected = (y_test / Y_tmp1_nor - hermval(x_tmp / popt[2], [0, 0, 0, array_hermite[3], 0, array_hermite[5]])) * Y_tmp1_nor
    popt2, pcov2 = curve_fit( gaussian, x_tmp, y_corrected,  sigma = Y_tmp1_err)
    array2[n]   = popt2[1]

plt.figure();
plt.plot(MJD, RV_HARPS-RVC, '.', MJD, array1, 'o', MJD, array2, '^')
plt.xlabel('MJD')
plt.ylabel('RV')
plt.title('Fitting up to 5th order')


'''
    
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