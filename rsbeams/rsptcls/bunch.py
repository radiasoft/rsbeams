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

Nonlinear waterbag and truncated Gaussian functions modified from C. Mitchell of LBNL.

Author: Nathan Cook
Date Created: 6/14/2018
Last Updated: 6/25/2018


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
        quiet (Boolean): boolean describing whether to use exact centroid injection, defaults to false.

    """
    
    def __init__(self, npart, dist = 'Gaussian', emitx = 1e-6, emity = 1e-6, betax=1., alphax = 0., 
                 betay =1., alphay=0., stdz=1., dpop=0, seed = None, quiet=False):
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
            quiet (Boolean): boolean describing whether to use exact centroid injection, defaults to false.
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
        
        #define quiet injection attribute
        self.quiet = quiet
          
        
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
                
                #Add 3 more particles if creating a quiet start
                if self.quiet:
                    self.exact_centroids(ptclCoords, phaseSpaceList)
                    ptclsMade += 3
        
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
        
    def exact_centroids(self, ptclCoords, phaseSpaceList):
        """Sets the centroid of the distribution in phase space to be zero - providing a quieter injection."""
        
        translation1 = np.array([-1,1,-1,1])
        translation2 = np.array([-1,-1,-1,-1])
        translation3 = np.array([1,-1,1,-1])

        ptcl1 = ptclCoords * translation1
        ptcl2 = ptclCoords * translation2
        ptcl3 = ptclCoords * translation3

        phaseSpaceList.append(ptcl1)
        phaseSpaceList.append(ptcl2)
        phaseSpaceList.append(ptcl3)
        
        
class NonlinearBunch(StandardBunch):
    
    """ 
    Derived class for generating matched distributions to a Danilov-Nagaitsev style nonlinear insert.
    
    Attributes:
        npart (int): the number of particles in the bunch
        dist (string): distribution type ['KV, 'waterbag', 'Gaussian'], defaults to KV
        emitx (float): RMS uncorrelated emittance in the x-px plane, defaults to 1 mm-mrad
        emity (float): RMS uncorrelated emittance in the y-py plane, defaults to 1 mm-mrad
        betax (float): the beta function where the bunch is being matched, defaults to 1
        alphax (float): one-half the derivative of the beta function, defaults to 0
        betay (float): the beta function where the bunch is being matched, defaults to 1
        alphay (float): one-half the derivative of the beta function, defaults to 0
        stdz (float): standard deviation in z-coordinate, defautls to 0
        dpop (float): standard deviation in delta-p/p0 coordinate, defaults to 0
        seed (float): seed for psuedorandom generator. Defaults to None, and is generate during initialization.
        quiet (Boolean): boolean describing whether to use exact centroid injection, defaults to false.
        t (float): nonlinear strength parameter for the insert. Defaults to 0.1 (unitless).
        c (float): the nonlinear aperture parameter (m^-1/2), defining poles in the x-axis. Defaults to 0.01.
        cutoff (float): cutoff parameter for the nonlinear Gaussian distributoin, defaults to 4.

    """
    
    def __init__(self, npart, dist = 'KV', emitx = 1e-6, emity = 1e-6, betax=1., alphax = 0., 
                 betay =1., alphay=0., stdz=1., dpop=0, seed = None, queit=False, t = 0.1, c = 0.01, cutoff = 4):
        """
        
        Args:
            npart (int): the number of particles in the bunch
            dist (string): distribution type ['KV, 'waterbag', 'Gaussian'], defaults to KV
            emitx (float): RMS uncorrelated emittance in the x-px plane, defaults to 1 mm-mrad
            emity (float): RMS uncorrelated emittance in the y-py plane, defaults to 1 mm-mrad
            betax (float): the beta function where the bunch is being matched, defaults to 1
            alphax (float): one-half the derivative of the beta function, defaults to 0
            betay (float): the beta function where the bunch is being matched, defaults to 1
            alphay (float): one-half the derivative of the beta function, defaults to 0
            stdz (float): standard deviation in z-coordinate, defautls to 0
            dpop (float): standard deviation in delta-p/p0 coordinate, defaults to 0
            seed (float): seed for psuedorandom generator. Defaults to None, and is generate during initialization.
            quiet (Boolean): boolean describing whether to use exact centroid injection, defaults to false.
            t (float): nonlinear strength parameter for the insert. Defaults to 0.1 (unitless).
            c (float): the nonlinear aperture parameter (m^-1/2), defining poles in the x-axis. Defaults to 0.01.
            cutoff (float): cutoff parameter for the nonlinear Gaussian distributoin, defaults to 4.
        """
        
        super(NonlinearBunch,self).__init__(npart, dist, emitx, emity, betax, alphax, betay, alphay, stdz, dpop, seed)
        
        self._t = t
        self._c = c
        self._cutoff = cutoff
        
        
    #Define c property which maintains a positive definite value
    @property
    def c(self):
        return self._c
    @c.setter
    def c(self, cval):
        if cval < 0:
            raise ValueError("Aperture parameter c must be larger than 0.")     
        self._c = cval
        
            
    #Define t property which maintains a positive definite value
    @property
    def t(self):
        return self._t
    @t.setter
    def t(self, tval):
        if tval < 0:
            raise ValueError("Nonlinear strength parameter t must be larger than 0.")     
        self._t = tval
        
    #Define cutoff property which maintains a positive definite value
    @property
    def cutoff(self):
        return self._cutoff
    @cutoff.setter
    def cutoff(self, cutoff):
        if cval < 0:
            raise ValueError("Cutoff parameter for distributions must be larger than 0.")     
        self._cutoff = cutoff
        
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
        
        #betax should be the same as betay
        self._betay = bet
        self._gammay = (1 + self._alphax**2) / self._betax
        
    
    @property
    def betay(self):
        return self._betay
    @betax.setter
    def betay(self, bet):
        if bet < 0:
            raise ValueError("Beta must be larger than 0.")
        
        self._betay = bet
        self._gammay = (1 + self._alphay**2) / self._betay
        
        #betax and betay should be the same
        self._betax = bet
        self._gammax = (1 + self._alphax**2) / self._betax
        
        
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
    
    
    def set_transverse_coordinates(self, emitx = None, emity = None, betax = None, alphax = None, cutoff = None):
        """Define the arrays describing the longitudinal coordinates z, dpop"""
        
        #The nonlinear insert requires equal beta functions, therefore we only permit setting betax
        if betax is not None:
            self.betax = betax
            
        if alphax is not None:
            self.alphax = alphax

        if emitx is not None:
            self.emitx = emitx
        if emity is not None:
            self.emity = emity
        if cutoff is not None:
            self.cutoff = cutoff
         
        if self.dist == 'KV':
            self.distribute_KV()
            
        elif self.dist == 'waterbag':
            self.distribute_waterbag()
            
        elif self.dist == 'Gaussian':
            self.distribute_Gaussian()
            
    def distribute_KV(self):
        """ 
        Generates a generalized KV distribution in 4D phase space using known bunch attributes. 
        Note that the KV distribution uniqely characterizes the bunch given a single emittance 
        and appropriate normalizing (Twiss) parameters. In the nonlinear case, there is a different
        relationship between the total emittance and the planar emittance.
        """
        
        assert (self.emitx == self.emity), "For a KV distribution, the planar emittances must be equal"
        
        #total emittance of the K-V distribution is equal to the planar emittance
        #this differs from the linear K-V distribution
        emit = self.emitx
        self.emit = emit
        
        # Generate some bounds on the transverse size to reduce waste in generating the bunch
        # Use the lemming method to find the maximum y
        y0 = np.sqrt(self.emit)
        
        yMax = newton(self.whatsleft, y0)
        
        #bounding the horizontal coordinate is difficult, but it should not exceed the pole
        xMax = self.c
        
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
                
                #Add 3 more particles if creating a quiet start
                if self.quiet:
                    self.exact_centroids(ptclCoords, phaseSpaceList)
                    ptclsMade += 3       
        
        self.particles[:,:4] = np.asarray(phaseSpaceList)
    
    
    def distribute_waterbag(self):
        """
        Generates a Waterbag distribution tailored to the elliptic potential. This version uses 
        a modified Jacobian factor dervied by Chad Mitchell to produce a smoothly varying distribution
        at the origin. A traditional Waterbag scheme produces an unphysical peak at the origin.
        
        The method of generating particle coordinates remains simular to that used for the K-V distribution.
        
        """
        # Generate particles by creating trials and finding particles with potential less than emittance, then assign the rest to momentum
        ptclsMade = 0
        phaseSpaceList = []
        while ptclsMade < self.npart:
            ranU = 0.0
            while ranU <= 0:
                ranU = random.random()
            
            # Generate some bounds on the transverse size to reduce waste in generating the bunch
            # Use the lemming method to find the maximum y
            trialH = np.sqrt(ranU)
            newH = self.emit*trialH
            y0 = np.sqrt(newH)
            #self.emittance = newH
            yMax = newton(self.whatsleft, y0) 
            
            #bounding the horizontal coordinate is difficult, but it should not exceed the pole
            xMax = self.c
            #xMax = yMax
            
            trialValue = 1e10
            while trialValue >= newH:
                xTrial = 2.*(0.5 - random.random())*xMax
                yTrial = 2.*(0.5 - random.random())*yMax
                trialValue = self.compute_potential(xTrial, yTrial)
            
            initialValue = trialValue
            if initialValue < newH:
                pMag = np.sqrt(2*(newH - initialValue))
                pDir = 2*np.pi* random.random()
                pxHat = pMag * np.cos(pDir)
                pyHat = pMag * np.sin(pDir)
                xReal = xTrial * np.sqrt(self.betax)
                yReal = yTrial * np.sqrt(self.betay)
                pxReal = (pxHat - self.alphax*xTrial)/np.sqrt(self.betax)
                pyReal = (pyHat - self.alphay*yTrial)/np.sqrt(self.betay)
                ptclCoords = np.array([xReal, pxReal, yReal, pyReal])
                phaseSpaceList.append(ptclCoords)
                ptclsMade += 1
                
                #Add 3 more particles if creating a quiet start
                if self.quiet:
                    self.exact_centroids(ptclCoords, phaseSpaceList)
                    ptclsMade += 3
            else:
                print "Initial value generated exceeds limiting H. Sampling new value."
        
        self.particles[:,:4] = np.asarray(phaseSpaceList)
    
    
    def distribute_Gaussian(self):
        """
        Generates a truncated Gaussian distribution in H-space for the elliptic potential according to:
            
            - P(H) = exp(-H/eps) for (0 < H/eps < L], 
            
            where H is the nonlinear Hamiltonian, eps is the emittance, and L is the cutoff parameter.
         
        The method of generating particle coordinates remains simular to that used for the K-V distribution.
        
        """
        # Copy the emittance value temporarily. It will be needed to reset the bunch attribute after fitting.
        bunch_emittance = self.emit
        
        # Generate particles by creating trials and finding particles with potential less than emittance, then assign the rest to momentum
        ptclsMade = 0
        phaseSpaceList = []
        while ptclsMade < self.npart:
            trialH = 1e10 #1.0e10
            while trialH > self.cutoff:        #Test against cutoff value
                ranU1 = 0.0
                ranU2 = 0.0
                while ranU1*ranU2 <= 0:
                    ranU1 = random.random()
                    ranU2 = random.random()
                trialH = -1.0*np.log(ranU1*ranU2)   #Generate an Erlang distribution in h
            
            
            # Generate some bounds on the transverse size to reduce waste in generating the bunch
            # Use the lemming method to find the maximum y            
            newH = bunch_emittance*trialH #must modify the original bunch emtitance here
            y0 = np.sqrt(newH)
            
            self.emit = newH #temporarily reset emittance for computing ymax on a particle by particle basis

            yMax = newton(self.whatsleft, y0) 
            
            #bounding the horizontal coordinate is difficult, but it should not exceed the pole
            xMax = self.c
            #xMax = yMax
            
            trialValue = 1e10
            while trialValue >= newH:
                xTrial = 2.*(0.5 - random.random())*xMax
                yTrial = 2.*(0.5 - random.random())*yMax
                trialValue = self.compute_potential(xTrial, yTrial)
            initialValue = trialValue
            if initialValue < newH:
                pMag = np.sqrt(2*(newH - initialValue))
                pDir = 2*np.pi* random.random()
                pxHat = pMag * np.cos(pDir)
                pyHat = pMag * np.sin(pDir)
                xReal = xTrial * np.sqrt(self.betax)
                yReal = yTrial * np.sqrt(self.betay)
                pxReal = (pxHat - self.alphax*xTrial)/np.sqrt(self.betax)
                pyReal = (pyHat - self.alphay*yTrial)/np.sqrt(self.betay)
                ptclCoords = np.array([xReal, pxReal, yReal, pyReal])
                phaseSpaceList.append(ptclCoords)
                ptclsMade += 1
                
                #Add 3 more particles if creating a quiet start
                if self.quiet:
                    self.exact_centroids(ptclCoords, phaseSpaceList)
                    ptclsMade += 3
            else:
                print "Initial value generated exceeds limiting H. Sampling new value."
        
        #Completed distribution, so reset emittance
        self.emit = bunch_emittance
        
        self.particles[:,:4] = np.asarray(phaseSpaceList)
    
    
    def compute_Hamiltonian(self, xHat, pxHat, yHat, pyHat):
        """Compute the Hamiltonian (1st invariant) for the integrable elliptic potential"""

        quadratic = 0.5 * (pxHat**2 + pyHat**2)

        hamiltonian = quadratic + self.compute_potential(xHat, yHat)
        
        return hamiltonian
        
    def compute_potential(self, xHat, yHat):
        """Compute the general potential for elliptic element with strength t"""
        
        quadratic = 0.5 * (xHat**2 + yHat**2)
        
        #compute default prefactors
        elliptic = 0.
        kfac = 1.
        
        #only recompute if t > 0
        if self._t != 0.:
            xN = xHat / self._c
            yN = yHat / self._c

            # Elliptic coordinates
            u = ( np.sqrt((xN + 1.)**2 + yN**2) +
                  np.sqrt((xN - 1.)**2 + yN**2) )/2.
            v = ( np.sqrt((xN + 1.)**2 + yN**2) -
                  np.sqrt((xN - 1.)**2 + yN**2) )/2.

            f2u = u * np.sqrt(u**2 - 1.) * np.arccosh(u)
            g2v = v * np.sqrt(1. - v**2) * (-np.pi/2 + np.arccos(v))

            kfac = -1.*self._t * self._c**2
            elliptic = (f2u + g2v) / (u**2 - v**2)

        potential = quadratic + kfac * elliptic
        
        return potential