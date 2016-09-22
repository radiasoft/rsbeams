"""
Phases.py is a module for performing analyses of phase advance and tunes for single particles or groups of particles.
As part of the rsbeams package, it should be code agnostic, and therefore can be imported and used for bunch diagnostics
for any code, so long as the particle coordinates adhere the analysis format.


Author: Nathan Cook
Date Created: 9/21/2016
Last Updated: 9/21/2016


"""

import sys
import os
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt


#Coords is a coordinate rubric for how particle phase space coordinates are assumed to be stored for particle beams
coords = {}
coords['x'] = 0
coords['xp'] = 1
coords['y'] = 2
coords['yp'] = 3
coords['cdt'] = 4
coords['dpop'] = 5
coords['id'] = 6


###################################### Phase Unwrap ###################################



###################################### FFT Frequency Analysis ###################################


def tune_fft(x_c, t_i = 1, t_f = None):
    '''
    Estimate the tune using an FFT of particle coordinates for N particles.
    Assuming we sample 1x per turn, then the tune resolution is 1/(t_f-t_i),
    the e.g. the reciprocal of the number of points in the sample interval.
    
    Arguments:
        x_c (ndArray): Nx1 array of particle coordinates in one plane
        t_i (Optional[int]): Index of x_c from which the fft begins. Defaults to 1.
        t_f (Optional[int]): Index of x_c at which the fft ends. Defaults to None.
    
    
    '''
    
    if not t_f:
        #if no end specified, then take interval to end of x_c
        t_f = len(x_c)
    
    
    num_used = len(x_c[t_s:t_f])
    tv = np.arange(num_used)*1.0/num_used
    sp = np.abs(np.fft.fft(x_c[t_s:t_f])) #take absolute value of every value of fft array

    smax = np.max(sp)
    m_ind = np.where(sp == smax)
    Q_guess =m_ind[0][0]*1./num_used
    if Q_guess > 0.5:
        Q_calc = 1.- Q_guess
    else:
        Q_calc = Q_guess
    return Q_calc