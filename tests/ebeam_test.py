from __future__ import absolute_import, division, print_function, unicode_literals
import pytest
import numpy
from rsbeams.physics import rsconst
from rsbeams.ptcl_beam import RsPtclBeam6D

def test_electron_beam():
    # specify physical properties of the beam
    num_ptcls = 100
    design_p_ev = 271
    total_charge_c = 3.05
    mass_ev = rsconst.m_e
    alpha_x = 1.3     # []
    beta_x = 21.1     # [m/rad]
    emit_x = 1.7e-06  # [m-rad]
    alpha_y = 1.3     # []
    beta_y = 21.1     # [m/rad]
    emit_y = 1.7e-06  # [m-rad]
    alpha_z = 1.3     # []
    beta_z = 21.1     # [m/rad]
    emit_z = 1.7e-06  # [m-rad]

    my_ebeam = RsPtclBeam6D.RsPtclBeam6D(num_ptcls, design_p_ev, \
                                         total_charge_c, mass_ev, \
                                         alpha_x, beta_x, emit_x, \
                                         alpha_y, beta_y, emit_y, \
                                         alpha_z, beta_z, emit_z )


test_electron_beam()
