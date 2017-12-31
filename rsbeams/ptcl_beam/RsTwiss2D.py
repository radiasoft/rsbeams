# -*- coding: utf-8 -*-
"""Data holding class for 2D Twiss parameters.

Original code taken from RadTrack project, https://github.com/radiasoft/radtrack
:copyright: Copyright (c) 2013 RadiaBeam Technologies, LLC. All Rights Reserved.

Subsequent mods are due to RadiaSoft,
:copyright: Copyright (c) 2017 Radiasoft LLC. All Rights Reserved.

:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""

class RsTwiss2D:
    """Data holding class for 2D Twiss parameters."""

    def __init__(self, alpha_rms, beta_rms, emit_rms):
        self.emit_rms  = emit_rms
        self.beta_rms  = beta_rms
        self.alpha_rms = alpha_rms
        self.gamma_rms = (1.0 + alpha_rms**2) / beta_rms

    def get_emit_rms(self):
        return self.emit_rms

    def get_beta_rms(self):
        return self.beta_rms

    def get_alpha_rms(self):
        return self.alpha_rms

    def get_gamma_rms(self):
        return self.gamma_rms
