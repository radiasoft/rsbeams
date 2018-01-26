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
from collections import OrderedDict

from rsbeams.rsptcls import RsDistrib6D
from rsbeams.rsptcls import RsTwiss2D
from rsbeams.rsphysics import rsconst

class RsPtclBeam6D:
    """Representation of a 6D charged particle distribution."""

    def __init__(self, num_ptcls, design_p_ev, total_charge_c, mass_ev, dist_type, max_rms_fac, \
                 alpha_x, beta_x, emit_x, alpha_y, beta_y, emit_y, alpha_z, beta_z, emit_z):

        # set the specified data
        self.design_p_ev = design_p_ev
        self.total_charge_c = total_charge_c
        self.mass_ev = mass_ev

        # instantiate local class variables
        self.twiss6d = OrderedDict( [('twiss_x', RsTwiss2D.RsTwiss2D(alpha_x, beta_x, emit_x)), \
                                     ('twiss_y', RsTwiss2D.RsTwiss2D(alpha_y, beta_y, emit_y)), \
                                     ('twiss_z', RsTwiss2D.RsTwiss2D(alpha_z, beta_z, emit_z))] )

        self.distrib6d = RsDistrib6D.RsDistrib6D(num_ptcls, dist_type, max_rms_fac)
        self.distrib6d.make_twiss_dist_6d(self.twiss6d, self.design_p_ev)

        return

    def get_design_p_ev(self):
        return self.design_p_ev

    def get_total_charge_c(self):
        return self.total_charge_c

    def get_mass_ev(self):
        return self.mass_ev

    def get_beta0_gamma0(self):
        return self.design_p_ev / self.mass_ev

    def get_gamma0(self):
        return math.sqrt(self.get_beta0_gamma0()**2 +1.)

    def get_beta0(self):
        return self.get_beta0_gamma0() / self.get_gamma0()

    # allowed values: 'twiss_x', 'twiss_y', 'twiss_z'
    def get_twiss2d_by_name(self, twiss_name):
        return self.twiss6d[twiss_name]

    def calc_twiss6d(self):
        self.distrib6d.calc_twiss6d(self.twiss6d)
        return

    def get_distrib6d(self):
        return self.distrib6d

    def get_peak_current_rms(self):
        rms_length = self.distrib6d.calc_rms_values_6d()[4]
        rms_time = rms_length / (self.get_beta0()*rsconst.c)
        return self.total_charge_c/rms_time
