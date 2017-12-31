from __future__ import absolute_import, division, print_function, unicode_literals
import pytest
import numpy
from rsbeams.ptcl_beam import RsDistrib6D

def test_unif_sphere():
    # initialization creates a 6 x N array of zeros
    my_n = 1000
    my_distrib = RsDistrib6D.RsDistrib6D(my_n)

    exception_thrown = False
    try:
        # type must be "uniform" or "gaussian"
        my_distrib.set_distrib_type('Uniform')
    except:
        exception_thrown = True

    assert(exception_thrown)

    my_distrib.set_distrib_type('uniform')
    my_distrib.make_unif_distrib()
#    print('\nmy_distrib =')
#    print(my_distrib.get_phase_space_6d().get_array_6d()[:,:])

    my_avg = my_distrib.calc_averages_6d()
    my_rms = my_distrib.calc_rms_values_6d()

#    print('\nmy_avg =', my_avg[:])
#    print('\nmy_rms =', my_rms[:])

    for i in range(6):
        assert(my_avg[i] < 0.04)
        assert(my_rms[i] < 0.48)
        assert(my_rms[i] > 0.42)

test_unif_sphere()
