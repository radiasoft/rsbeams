"""
nonlinear.py is a module containing tools for creating nonlinear lattice elements and corresponding components, 
specifically emphasizing the nonlinear integrable insert for use in IOTA and other integrable systems using 
the elliptic scheme defined by: https://journals.aps.org/prab/abstract/10.1103/PhysRevSTAB.13.084002
V. Danilov and S. Nagaitsev. Phys. Rev. ST Accel. Beams 13, 084002 (2010
Author: Nathan Cook
Date Created: 6/21/2018
Last Updated: 6/21/2018
"""

import numpy as np


class NonlinearInsert(object):
    """
    Class for generating and manipulating a nonlinear insert for use in the nonlinear integrable optics.
    Following the presciptions of the elliptic scheme defined by: 
       -V. Danilov and S. Nagaitsev. Phys. Rev. ST Accel. Beams 13, 084002 (2010)
       -https://journals.aps.org/prab/abstract/10.1103/PhysRevSTAB.13.084002    
    
    Instantiation of the nonlinear insert requires element requires specification of the length, phase, 
    aperture, strength, and number of slices comprising the insert. The required Twiss functions can be
    computed from this information, and hence insertion of the element into a lattice can be done by 
    comparing the requirements to the Twiss parameters obtained from the lattice.
    
    Attributes:
        length (float): the length of the nonlinear insert in meters
        phase (float): the phase advance modulo 2pi through the nonlinear insert
        t (float): the dimensionless nonlinear strength parameter. Defaults to 0.1
        c (float): the nonlinear aperture parameter (m^-1/2), defining poles in the x-axis. Defaults to 0.01.
        num_slices (int): the number of piecewise constant segements used to construct the insert
        
        s_vals (ndArray): array of relative s-values providing the center of each segment of the nonlinear insert
        knll (ndArray): array of gradient values for each nonlinear segment element (nllens)
        cnll (ndArray): array of aperture parameters for each nonlinear segment element (nllens)
       
    """
    
    def __init__(self, length, phase, t = 0.1, c = 0.01, num_slices = 20):
    
        """
        Arguments:
            length (float): the length of the nonlinear insert in meters
            phase (float): the phase advance modulo 2pi through the nonlinear insert
            t (float): the dimensionless nonlinear strength parameter. Defaults to 0.1
            c (float): the nonlinear aperture parameter (m^-1/2), defining poles in the x-axis. Defaults to 0.01.
            num_slices (int): the number of piecewise constant segements used to construct the insert
        """
        self.length = length
        self.phase = phase
        self.t = t
        self._c = c
        self.num_slices = num_slices
    
        
    #Define c property which maintains a positive definite value
    @property
    def c(self):
        return self._c
    @c.setter
    def c(self, cval):
        if cval < 0:
            raise ValueError("Aperture parameter c must be larger than 0.")     
        self._c = c
        
    
    def generate_sequence(self):
        """Generates arrays containing the knll and cnll values for each nllens element"""
        
        #Define the focal length of the insert using the phase advance and length
        f0 = self.length/4.0*(1.0+1.0/np.tan(np.pi*self.phase)**2)
        
        #define array of s-values
        start = (self.length/self.num_slices)*0.5
        end = self.length - start
        #Make an attribute as they could be useful for constructing the mad-x sequence
        self.s_vals = np.linspace(start,end,self.num_slices) 
        
        #set the initial beta value to help compare to lattice functions
        self.beta0 = self.length*(1.-0.0*(self.length)/self.length/f0)/np.sqrt(1.0-(1.0-self.length/2.0/f0)**2)
        
        #set the beta functions as an attribute for comparing against other lattice functions
        bn = self.length*(1-self.s_vals*(self.length-self.s_vals)/self.length/f0)/np.sqrt(1.0-(1.0-self.length/2.0/f0)**2)
        self.betas = bn
        
        knn = self.t*self.length/self.num_slices/bn**2
        cnll = self.c*np.sqrt(bn)
        knll = knn*cnll**2
        
        #Now set the knll and cnll parameters for each nllens object
        self.knll = knll
        self.cnll = cnll
        
    def validate_sequence(self, beta_values):
        """
        Checks the predicted beta functions for the specified sequence against known lattice functions.
        This function is not currently implemented.
        
        Arguments:
            beta_values: an array of beta x/y values for the underlying "bare" lattice
            
        Returns:
            check: a boolean signifying that the lattice correctly permits the sequence defined in the lattice
        """
        
        return
        
            
    def create_madx(self):
        """
        Return a sequence of madx elements representing the insert, represented as strings:
        
            ["elem 1 string", "elem 2 string", ...]
        
        Returns:
            MADX_elements: list of strings describing nllens elements needed to construct inset
        """
        
        MADX_elements = []
        
        for ind in range(len(self.knll)):
            #Loop through element and construct nllens with strengths and apertures
            elem = "nllens, knll = {}, cnll = {};".format(self.knll[ind], self.cnll[ind])
            MADX_elements.append(elem)
            
            return MADX_elements
