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

def phase_advance(p1, p2, clockwise=True):
    '''
    Returns the angle between two vectors.
    
    This script is used to compute an effective phase advance by comparing an initial and final position in phase space.
    
    Arguments:
        p1 (ndarray): Array containing phase space coordinates (e.g. [x, x'])
        p2 (ndarray): Array containing phase space coordinates (e.g. [x, x'])
        clockwise (optional[Bool]): If true, assumes rotation is clockwise. Defaults to true.
        
    Returns:
        guess_angle (float): The angle between the input vectors, given in radians. Assumes clockwise rotation.
    
    '''
    #use the quadrant-dependent arcantangent to get the effective angle (measured clockwise) for each vector
    #and take the difference of the two.
    guess  = np.arctan2(p1[1],p1[0]) - np.arctan2(p2[1],p2[0])
    
    if guess < 0:
        advance = 2*np.pi + guess
    else:
        advance = guess
        
    if not clockwise:
        advance = 2*np.pi - advance
        
    return advance

###################################### FFT Frequency Analysis ###################################


def tune_fft(x_c, t_i=1, t_f=None):
    """
    Estimate the tune using an FFT of particle coordinates for N particles.
    Assuming we sample 1x per turn, then the tune resolution is 1/(t_f-t_i),
    the e.g. the reciprocal of the number of points in the sample interval.
    
    Arguments:
        x_c (ndArray): Nx1 array of particle coordinates in one plane
        t_i (Optional[int]): Index of x_c from which the fft begins. Defaults to 1.
        t_f (Optional[int]): Index of x_c at which the fft ends. Defaults to None.
    """
    
    if not t_f:
        # if no end specified, then take interval to end of x_c
        t_f = len(x_c)

    num_used = len(x_c[t_i:t_f])
    sp = np.abs(np.fft.fft(x_c[t_i:t_f]))  # take absolute value of every value of fft array

    m_ind = np.argmax(sp)
    Q_guess =m_ind * 1. / num_used
    if Q_guess > 0.5:
        Q_calc = 1. - Q_guess
    else:
        Q_calc = Q_guess
    return Q_calc
