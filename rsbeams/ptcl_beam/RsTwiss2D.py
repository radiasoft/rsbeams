# -*- coding: utf-8 -*-
"""Data holding class for 2D Twiss parameters.

Original code taken from RadTrack project, https://github.com/radiasoft/radtrack
:copyright: Copyright (c) 2013 RadiaBeam Technologies, LLC. All Rights Reserved.

Subsequent mods are due to RadiaSoft,
:copyright: Copyright (c) 2017 Radiasoft LLC. All Rights Reserved.

:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""

class RbTwiss2D:
    """Data holding class for 2D Twiss parameters."""

    def __init__(self, alphaRMS, betaRMS, emitRMS):
        self.emitRMS  = emitRMS
        self.betaRMS  = betaRMS
        self.alphaRMS = alphaRMS
        self.gammaRMS = (1.0 + alphaRMS**2) / betaRMS

    def getEmitRMS(self):
        return self.emitRMS

    def getBetaRMS(self):
        return self.betaRMS

    def getAlphaRMS(self):
        return self.alphaRMS

    def getGammaRMS(self):
        return self.gammaRMS
