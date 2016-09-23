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
    Returns the phase advance between two coordinates (2-vectors).
    
    Computes an effective phase advance by comparing an initial and final position in phase space.
    
    Arguments:
        p1 (ndarray): Array containing phase space coordinates (e.g. [x, x'])
        p2 (ndarray): Array containing phase space coordinates (e.g. [x, x'])
        clockwise (optional[Bool]): If true, assumes rotation is clockwise. Defaults to true.
        
    Returns:
        guess_angle (float): The angle between the input vectors, given in radians. Assumes clockwise rotation.
    
    '''
    
    norm1 = np.linalg.norm(p1)
    norm2 = np.linalg.norm(p2)
    
    #Determine the quadrant - the integer corresponds to the quadrant
    q1 = 1*((p1[0] > 0) & (p1[1] >= 0)) + 2*((p1[0] <= 0) & (p1[1] > 0)) + 3*((p1[0] < 0) & (p1[1] <= 0)) + 4*((p1[0] >= 0) & (p1[1] < 0))
    q2 = 1*((p2[0] > 0) & (p2[1] >= 0)) + 2*((p2[0] <= 0) & (p2[1] > 0)) + 3*((p2[0] < 0) & (p2[1] <= 0)) + 4*((p2[0] >= 0) & (p2[1] < 0))

    
    #calculate dot product
    product = np.dot(p1,p2)
    guess_angle = np.arccos(product/(norm1*norm2))
    
    #determine if phase advance > 180 degrees
    #For the 3,4 -> 1,2 cases, always want the negative x value to have larger magnitude
    if (q2-q1 >= 2):
        #rotation > 180
        greater = True
    elif (q1==3) and (q2==1) and (np.abs(p1[0])>p2[0]):
        #rotation > 180 because x1 > x2 for q3 -> q1
        greater = True
    elif (q1==4) and (q2==2) and (np.abs(p2[0])>p1[0]):
        #rotation > 180 because x2 > x1 for q4 -> q2
        greater = True
    else:
        greater = False
    
    #if phase advance > 180 for counterclockwise rotation, return 2.*np.pi - guess_angle
    if greater:
        if clockwise:
            return guess_angle
        else:
            return 2*np.pi- guess_angle    
    else:
        if clockwise:
            return 2*np.pi- guess_angle
        else:
            return guess_angle

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
