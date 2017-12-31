# -*- coding: utf-8 -*-
"""Representation of a 6D charged particle distribution.

We follow the conventions of the Elegant code from Argonne National Lab:
  The fundamental data structure is N rows by 6 columns
  N is the number of macroparticles
  We assume 6 phase space coordinates/momenta
     x [m], xp=px/p0 [rad], y [m], yp=py/p0 [rad], s [m], (p-p0)/p0 [rad]
  Here, s is the total distance traveled along the accelerator axis,
        and p0 is the design momentum [eV/c]

Original code taken from RadTrack project, https://github.com/radiasoft/radtrack
:copyright: Copyright (c) 2014 RadiaBeam Technologies, LLC. All Rights Reserved.

Subsequent mods are due to RadiaSoft,
:copyright: Copyright (c) 2017 Radiasoft LLC. All Rights Reserved.

:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""
import math
import scipy.constants
from collections import OrderedDict

from rsbeams.ptcl_beam import RsDistrib6D
from rsbeams.ptcl_beam import RsTwiss2D

class RsPtclBeam6D:
"""Representation of a 6D charged particle distribution."""

    def __init__(self, numPtcls, designMomentumEV, totalCharge, massEV,
                 alphaX, betaX, emitX, alphaY, betaY, emitY, alphaZ, betaZ, emitZ):

        # set the specified data
        self.designMomentumEV = designMomentumEV
        self.totalCharge = totalCharge
        self.massEV = massEV

        # instantiate local class variables
        self.distrib6D = RsDistrib6D.RsDistrib6D(numPtcls)
        self.twissParams6D = OrderedDict( [('twissX', RsTwiss2D(alphaX, betaX, emitX)), \
                                           ('twissY', RsTwiss2D(alphaX, betaX, emitX)), \
                                           ('twissZ', RsTwiss2D(alphaX, betaX, emitX))] )

        # define useful temporary variables
        self.rt2opi = math.sqrt(2./math.pi)

        # specify physical constants
        self.c     = scipy.constants.c          # speed of light [m/s]
        self.cSq   = self.c**2            # speed of light squared
        self.cInv  = 1./self.c            # one over the speed of light
        self.mu0   = scipy.constants.mu_0    # permeability of free space
        self.eps0  = scipy.constants.epsilon_0 # permittivity of free space
        self.eMass   = scipy.constants.m_e     # electron mass [kG]
        self.eCharge = scipy.constants.e   # elementary charge [C]
        self.eMassEV = self.eMass*self.cSq/self.eCharge  # eMass [eV]

        return

    def getDesignMomentumEV(self):
        return self.designMomentumEV

    def setDesignMomentumEV(self, designMomentumEV):
        if (designMomentumEV > 0.):
            self.designMomentumEV = designMomentumEV
        else:
            message = 'ERROR!  designMomentumEV <= 0.: ' + str(designMomentumEV)
            raise Exception(message)
        return

    def getTotalCharge(self):
        return self.totalCharge

    def setTotalCharge(self, totalCharge):
        if (totalCharge > 0.):
            self.totalCharge = totalCharge
        else:
            message = 'ERROR!  totalCharge <= 0.: ' + str(totalCharge)
            raise Exception(message)
        return

    def getMassEV(self):
        return self.massEV

    def setMassEV(self, massEV):
        if (massEV > 0.):
            self.massEV = massEV
        else:
            message = 'ERROR!  massEV <= 0.: ' + str(massEV)
            raise Exception(message)
        return

    def getBeta0Gamma0(self):
        return self.designMomentumEV / self.massEV

    def getGamma0(self):
        return math.sqrt(self.getBeta0Gamma0()**2 +1.)

    def getBeta0(self):
        return self.getBeta0Gamma0() / self.getGamma0()

    # allowed values: 'twissX', 'twissY', 'twissZ'
    def getTwissParamsByName2D(self, twissName):
        return self.twissParams6D[twissName]

    def setTwissParamsByName2D(self, alpha, beta, emittance, twissName):
        self.twissParams6D[twissName] = RsTwiss2D.RsTwiss2D(alpha, beta, emittance)
        return

    def calcTwissParams6D(self):
        self.distrib6D.calcTwissParams6D(self.twissParams6D)
        return

    def getDistrib6D(self):
        return self.distrib6D

    def makePtclPhaseSpace6D(self, meanMomentum):
        self.distrib6D.makeTwissDist6D(self.twissParams6D, meanMomentum)
        return

    def getCurrent(self):
        s=self.distrib6D.calcRmsValues6D()[4]
        t=s/(self.getBeta0()*self.c)
        I = self.totalCharge/t
        return I
