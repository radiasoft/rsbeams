from __future__ import absolute_import, division, print_function, unicode_literals
import pytest
import numpy
from rsbeams.physics import rsconst
from rsbeams.statistics import stats6d
from rsbeams.ptcl_beam import RsPtclBeam6D

def test_beam_gen():

    # specify physical properties of the beam
    num_ptcls = 1000
    design_p_ev = 271
    total_charge_c = 3.05
    mass_ev = rsconst.m_e

    dist_type = 'gaussian'
    max_rms_fac = 4.9

    alpha_x = 1.3     # []
    beta_x = 21.1     # [m/rad]
    emit_x = 1.7e-06  # [m-rad]

    alpha_y = -3.01   # []
    beta_y = 9.07     # [m/rad]
    emit_y = 0.44e-6  # [m-rad]

    alpha_z = 0       # []
    beta_z =  1.55    # [m/rad]
    emit_z = 3.1e-05  # [m-rad]

    my_ebeam = RsPtclBeam6D.RsPtclBeam6D(num_ptcls, design_p_ev, \
                                         total_charge_c, mass_ev, \
                                         dist_type, max_rms_fac, \
                                         alpha_x, beta_x, emit_x, \
                                         alpha_y, beta_y, emit_y, \
                                         alpha_z, beta_z, emit_z )

    beta_gamma = my_ebeam.get_beta0_gamma0()
    gamma0 = my_ebeam.get_gamma0()
    beta0 = my_ebeam.get_beta0()

    print('\nbeta0, gamma0 = ', beta0, ', ', gamma0)

    my_twiss_x = my_ebeam.get_twiss2d_by_name('twiss_x')
    my_twiss_y = my_ebeam.get_twiss2d_by_name('twiss_y')
    my_twiss_z = my_ebeam.get_twiss2d_by_name('twiss_z')

#    print('alpha_x = ', my_twiss_x.get_alpha_rms())
#    print('beta_x = ', my_twiss_x.get_beta_rms())
#    print('emit_x = ', my_twiss_x.get_emit_rms())

    assert(my_twiss_x.get_alpha_rms() == alpha_x)
    assert(my_twiss_x.get_beta_rms() == beta_x)
    assert(my_twiss_x.get_emit_rms() == emit_x)

#    print('alpha_y = ', my_twiss_y.get_alpha_rms())
#    print('beta_y = ', my_twiss_y.get_beta_rms())
#    print('emit_y = ', my_twiss_y.get_emit_rms())

    assert(my_twiss_y.get_alpha_rms() == alpha_y)
    assert(my_twiss_y.get_beta_rms() == beta_y)
    assert(my_twiss_y.get_emit_rms() == emit_y)

    print('alpha_z = ', my_twiss_z.get_alpha_rms())
    print('beta_z = ', my_twiss_z.get_beta_rms())
    print('emit_z = ', my_twiss_z.get_emit_rms())

    assert(my_twiss_z.get_alpha_rms() == alpha_z)
    assert(my_twiss_z.get_beta_rms() == beta_z)
    assert(my_twiss_z.get_emit_rms() == emit_z)

    my_ebeam.calc_twiss6d()

    new_twiss_x = my_ebeam.get_twiss2d_by_name('twiss_x')
    new_twiss_y = my_ebeam.get_twiss2d_by_name('twiss_y')
    new_twiss_z = my_ebeam.get_twiss2d_by_name('twiss_z')

#    print('alpha_x = ', new_twiss_x.get_alpha_rms())
#    print('beta_x = ', new_twiss_x.get_beta_rms())
#    print('emit_x = ', new_twiss_x.get_emit_rms())

    assert(stats6d.specify_significant_figures(new_twiss_x.get_alpha_rms(),3) == \
           stats6d.specify_significant_figures(alpha_x,3))
    assert(stats6d.specify_significant_figures(new_twiss_x.get_beta_rms(),3) == \
           stats6d.specify_significant_figures(beta_x,3))
    assert(stats6d.specify_significant_figures(new_twiss_x.get_emit_rms(),3) == \
           stats6d.specify_significant_figures(emit_x,3))

#    print('alpha_y = ', new_twiss_y.get_alpha_rms())
#    print('beta_y = ', new_twiss_y.get_beta_rms())
#    print('emit_y = ', new_twiss_y.get_emit_rms())

    assert(stats6d.specify_significant_figures(new_twiss_y.get_alpha_rms(),3) == \
           stats6d.specify_significant_figures(alpha_y,3))
    assert(stats6d.specify_significant_figures(new_twiss_y.get_beta_rms(),3) == \
           stats6d.specify_significant_figures(beta_y,3))
    assert(stats6d.specify_significant_figures(new_twiss_y.get_emit_rms(),3) == \
           stats6d.specify_significant_figures(emit_y,3))

    print('alpha_z = ', new_twiss_z.get_alpha_rms())
    print('beta_z = ', new_twiss_z.get_beta_rms())
    print('emit_z = ', new_twiss_z.get_emit_rms())

#    assert(stats6d.specify_significant_figures(new_twiss_z.get_alpha_rms(),3) == \
#           stats6d.specify_significant_figures(alpha_z,3))
#    assert(stats6d.specify_significant_figures(new_twiss_z.get_beta_rms(),3) == \
#           stats6d.specify_significant_figures(beta_z,3))
#    assert(stats6d.specify_significant_figures(new_twiss_z.get_emit_rms(),3) == \
#           stats6d.specify_significant_figures(emit_z,3))

test_beam_gen()
