#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 16:57:03 2017

@author: jzhao
"""

"""
    Extract CCF: .fits -> .dat
"""

###########
# Updates #
###########
# bc = hdulist[0].header['HIERARCH ESO DRS BERV'] @26/07/17

#############################################

import sys
from astropy.io import fits
import glob
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
from math import sqrt
import math

############# DEFINE FUNCTIONS ##############
def gaussian(x, a, mu, sigma, C):
    val = a / ( sqrt(2*math.pi) * sigma ) * np.exp(-(x - mu)**2 / sigma**2) + C
    return val    

# Root mean square
def rms(num):
    return sqrt(sum(n*n for n in num)/len(num))
#############################################

FILE        = glob.glob('../ccf_fits/*fits')
N           = len(FILE)
N_start     = 0
N_end       = N
n_file      = N_end - N_start
MJD         = np.zeros(n_file)
RV_g        = np.zeros(n_file)

x           = np.arange(-10, -1+0.1, 0.1)                                # over sampling to 0.1 km/s [-10.2, -0.8]
y           = np.zeros(len(x))

plt.figure()

for n in range(N_start, N_end):
    
    # progress bar #
    sys.stdout.write('\r')
    sys.stdout.write("[%-50s] %d%%" % ('='*int((n+1-N_start)*50./(N_end-N_start)), int((n+1-N_start)*100./(N_end-N_start))))
    sys.stdout.flush()    

    
    hdulist     = fits.open(FILE[n])
    bc          = hdulist[0].header['HIERARCH ESO DRS BERV'] / 1000             # barycentric RV, originally in m/s
    v0          = hdulist[0].header['CRVAL1'] - bc                              # velocity on the left (N_starting point)
    # v0          = hdulist[0].header['CRVAL1']  
    MJD[n]      = hdulist[0].header['MJD-OBS']
    
    CCF         = hdulist[0].data                                               # ccf 2-d array
    ccf         = CCF[- 1, :]                                                   # ccf 1-d array (whole range)
    delta_v     = hdulist[0].header['CDELT1']                                   # velocity grid size 
    v           = v0 + np.arange(CCF.shape[1]) * delta_v                        # velocity array (whole range)

    if n == 0:                                                                  # not used
        idx_min     = ccf == min(ccf)
        v_min       = v[idx_min]

    f           = CubicSpline(v, ccf / ccf.max())
    y           = f(x)
    popt, pcov  = curve_fit( gaussian, x, y, [-1.76, -5.467, 2.5, 1])
#    print((popt[1] - v_min)*1000)
    RV_g[n]     = popt[1]
    x_new       = x
    y_new       = (y - popt[3]) / popt[0]
    # plt.plot(x_new, y_new, '-')
    
    writefile   = '../ccf_dat/ccf' + str(n) + '.dat'
    np.savetxt(writefile, y_new)

plt.plot(MJD, (RV_g - np.mean(RV_g))*1000, '.')    
plt.show()