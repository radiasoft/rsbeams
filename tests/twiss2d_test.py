from __future__ import absolute_import, division, print_function, unicode_literals
import pytest
import numpy
from rsbeams.rsstats import stats6d
from rsbeams.rsptcls import RsTwiss2D

def test_twiss2d():
    my_alpha = -0.1
    my_beta  = 13.
    my_emit  = 1.17e-06
    my_twiss = RsTwiss2D.RsTwiss2D(my_alpha, my_beta, my_emit)
    my_gamma = (1.0 + my_alpha**2) / my_beta
    test_gam = my_twiss.get_gamma_rms()

#    print('my_gamma = ', my_gamma)
#    print('test_gam = ', test_gam)

    assert stats6d.specify_significant_figures(my_gamma, 4) == \
           stats6d.specify_significant_figures(test_gam, 4)

# test_twiss2d()
