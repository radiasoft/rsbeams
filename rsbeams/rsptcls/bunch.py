"""
bunch.py is a module for computing bunch distributions for use with different tracking codes. The StandardBunch() class
defines beam phase space coordinates and permits fitting of different distributions for a traditional linear system
with given Twiss parameters. Currently, the class supports unit conventions for the Synergia tracking code, as
discussed in the Synergia 2.1 documentation: http://compacc.fnal.gov/~amundson/html/units.html

Unit Conventions:
    -x,y    : m
    -xp,yp  : unitless (px/ptotal)
    -z      : m (c*dt)
    -zp     : unitless (dpz/ptotal)

A derived class will soon be added for matching to a nonlinear integrable insert for use in IOTA and other integrable
systems using the elliptic scheme defined by: https://journals.aps.org/prab/abstract/10.1103/PhysRevSTAB.13.084002
V. Danilov and S. Nagaitsev. Phys. Rev. ST Accel. Beams 13, 084002 (2010


Author: Nathan Cook
Date Created: 6/14/2018
Last Updated: 6/14/2018


"""


import numpy as np
import random
from scipy.optimize import newton

class StandardBunch(object):
    
    """ 
    Generic class for generating traditional bunch distributions. 
    Generates a numpy array for easy output/input into other codes.
    
    Attributes:
        npart (int): the number of particles in the bunch
        dist (string): the distribution type being used, defaults to Gaussian
        emitx (float): RMS uncorrelated emittance in the x-px plane, defaults to 1 mm-mrad
        emity (float): RMS uncorrelated emittance in the y-py plane, defaults to 1 mm-mrad
        betax (float): the beta function where the bunch is being matched, defaults to 1
        alphax (float): one-half the derivative of the beta function, defaults to 0
        betay (float): the beta function where the bunch is being matched, defaults to 1
        alphay (float): one-half the derivative of the beta function, defaults to 0
        stdz (float): standard deviation in z-coordinate, defautls to 0
        dpop (float): standard deviation in delta-p/p0 coordinate, defaults to 0
        seed (float): seed for psuedorandom generator. Defaults to None, and is generate during initialization.

    """
    
    def __init__(self, npart, dist = 'Gaussian', emitx = 1e-6, emity = 1e-6, betax=1., alphax = 0., 
                 betay =1., alphay=0., stdz=1., dpop=0, seed = None):
        """
        
        Args:
            npart (int): the number of particles in the bunch
            dist (string): the distribution type being used, defaults to Gaussian
            emitx (float): RMS uncorrelated emittance in the x-px plane, defaults to 1 mm-mrad
            emity (float): RMS uncorrelated emittance in the y-py plane, defaults to 1 mm-mrad
            betax (float): the beta function where the bunch is being matched, defaults to 1
            alphax (float): one-half the derivative of the beta function, defaults to 0
            betay (float): the beta function where the bunch is being matched, defaults to 1
            alphay (float): one-half the derivative of the beta function, defaults to 0
            stdz (float): standard deviation in z-coordinate, defautls to 0
            dpop (float): standard deviation in delta-p/p0 coordinate, defaults to 0
            seed (float): seed for psuedorandom generator. Defaults to None, and is generate during initialization.
        """
        
        self.npart = npart
        self.dist = dist
        
        #set emittance parameters in each plane
        self.emitx = emitx
        self.emity = emity
        #total emittance - this varies with distribution but is useful to set a default 
        #it will be adjusted for specific distributions
        self.emit = emitx + emity
        
        #set Twiss parameters as "private" attributes because setters will be redefined to avoid conflicts
        self._betax  = betax
        self._alphax = alphax
        self._betay = betay
        self._alphay = alphay
        
        #use beta/alpha values to compute initial gamma
        self._gammax = (1 + alphax**2) / betax
        self._gammay = (1 + alphay**2) / betay
 

        self.stdz = stdz
        self.dpop = dpop
        
        #define seed
        if seed is not None:
            self.seed = seed
        else:
            self.seed = random.seed()
        
        #create particles array
        self.particles = np.zeros((npart,7))
        
        #define particle IDs
        self.particles[:,6] = np.arange(npart)
        
          
        
    #Define beta and alpha properties which automatically update gamma
    @property
    def betax(self):
        return self._betax
    @betax.setter
    def betax(self, bet):
        if bet < 0:
            raise ValueError("Beta must be larger than 0.")
        
        self._betax = bet
        self._gammax = (1 + self._alphax**2) / self._betax
    
    @property
    def betay(self):
        return self._betay
    @betax.setter
    def betay(self, bet):
        if bet < 0:
            raise ValueError("Beta must be larger than 0.")
        
        self._betay = bet
        self._gammay = (1 + self._alphay**2) / self._betay
        
    @property
    def alphax(self):
        return self._alphax
    @alphax.setter
    def alphax(self, alph):
        
        self._alphax = alph
        self._gammax = (1 + self._alphax**2) / self._betax
        
    @property
    def alphay(self):
        return self._alphay
    @alphay.setter
    def alphay(self, alph):
        
        self._alphay = alph
        self._gammay = (1 + self._alphay**2) / self._betay
    
    
    def set_longitudinal_coordinates(self, stdz = None, dpop = None):
        """Define the arrays describing the longitudinal coordinates z, dpop"""
        
        if stdz is not None:
            self.stdz = stdz
        if dpop is not None:
            self.dpop = dpop
        
        self.particles[:,4] = np.random.randn(self.npart)*self.stdz #set z coordinate
        self.particles[:,5] = np.random.randn(self.npart)*self.dpop #set dpop coordinate 
    
    def set_transverse_coordinates(self, emitx = None, emity = None, betax = None, alphax = None, betay = None, alphay = None):
        """Define the arrays describing the longitudinal coordinates z, dpop"""
        
        if betax is not None:
            self.betax = betax
        if alphax is not None:
            self.alphax = alphax
        if betay is not None:
            self.betax = betay        
        if alphay is not None:
            self.alphay = alphay
        if emitx is not None:
            self.emitx = emitx
        if emity is not None:
            self.emity = emity
        
        if self.dist == 'Gaussian':
            self.distribute_Gaussian()
            
        elif self.dist == 'KV':
            self.distribute_KV()
    
    
    def distribute_Gaussian(self):
        """ Generates an uncorrelated Gaussian distribution in 4D phase space using known bunch attributes"""
        
        sigma_x = np.sqrt(self.emitx*self._betax)
        sigma_xp = np.sqrt(self.emitx*self._gammax)
        
        sigma_y = np.sqrt(self.emity*self._betay)
        sigma_yp = np.sqrt(self.emity*self._gammay)
        
        self.particles[:,0] = np.random.randn(self.npart)*sigma_x #set x-coordinates
        self.particles[:,1] = np.random.randn(self.npart)*sigma_xp #set xp-coordinates
        self.particles[:,2] = np.random.randn(self.npart)*sigma_y #set y-coordinates
        self.particles[:,3] = np.random.randn(self.npart)*sigma_yp #set yp-coordinates
        
        
    def distribute_KV(self):
        """ 
        Generate a KV distribution in 4D phase space using known bunch attributes. Note that the KV distribution
        uniqely characterizes the bunch given a single emittance and appropriate normalizing (Twiss) parameters.
        """
        
        assert (self.emitx == self.emity), "For a KV distribution, the planar emittances must be equal"
        
        #total emittance of the K-V distribution is 4 times the planar emittance
        emit = 4.*self.emitx
        self.emit = emit
        
        # Generate some bounds on the transverse size to reduce waste in generating the bunch
        # Use the lemming method to find the maximum y
        y0 = np.sqrt(self.emit)
        
        yMax = newton(self.whatsleft, y0)        
        xMax = yMax
        
        # Generate particles by creating trials and finding particles with potential less than emittance, 
        # then assign the rest to momentum
        ptclsMade = 0
        phaseSpaceList = []
        
        while ptclsMade < self.npart:
            #Note that the particle coordinates here are distributed in normal coordinates
            xTrial = 2.*(0.5 - random.random())*xMax
            yTrial = 2.*(0.5 - random.random())*yMax
            trialValue = self.compute_potential(xTrial, yTrial)
            if trialValue < self.emit:
                
                pMag = np.sqrt(2.*(self.emit - trialValue))
                pDir = 2.*np.pi * random.random()
                pxHat = pMag * np.cos(pDir)
                pyHat = pMag * np.sin(pDir)
                
                xReal = xTrial * np.sqrt(self.betax)
                yReal = yTrial * np.sqrt(self.betay)
                
                #We want to provide the user with standard (non-normal) coordinates
                pxReal = (pxHat - self.alphax*xTrial)/np.sqrt(self.betax)
                pyReal = (pyHat - self.alphay*yTrial)/np.sqrt(self.betay)
                
                ptclCoords = np.array([xReal, pxReal, yReal, pyReal])
                phaseSpaceList.append(ptclCoords)
                ptclsMade += 1        
        
        self.particles[:,:4] = np.asarray(phaseSpaceList)
    
    
    
    
    def print_Twiss(self):
        """Print the Twiss parameters for the lattice being used to compute the bunch coordinates"""
        
        print "Twiss parameters in use:"
        print "betax : {}".format(self._betax)
        print "betay : {}".format(self._betay)
        print "alphax : {}".format(self._alphax)
        print "alphay : {}".format(self._alphay)        
        print "gammax : {}".format(self._gammax)
        print "gammay : {}".format(self._gammay)
        
    
    def compute_Hamiltonian(self, xHat, pxHat, yHat, pyHat):
        """Compute the Hamiltonian (C-S invariant) for the potential"""
        
        hamiltonian = 0.5*(pxHat**2 + pyHat**2) + 0.5 *(xHat**2 + yHat**2)
        
        return hamiltonian
        
    def compute_potential(self, xHat, yHat):
        """Compute the general potential"""
        potential = 0.5*(xHat**2 + yHat**2)
        
        return potential
        
    def whatsleft(self, yHat):
        """Return the difference btween the emittance and potential"""
        
        return self.emit - self.compute_potential(0, yHat)