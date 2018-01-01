from __future__ import absolute_import, division, print_function, unicode_literals
import pytest
import numpy
from rsbeams.statistics import stats6d
from rsbeams.ptcl_beam import RsDistrib6D

def test_unif_sphere():
    # initialization creates a 6 x N array of zeros
    my_n = 1000
    my_type = 'uniform'
    my_rms_fac = 4.7
    my_distrib = RsDistrib6D.RsDistrib6D(my_n, my_type, my_rms_fac)

#    print('\nmy_distrib =')
#    print(my_distrib.get_phase_space_6d().get_array_6d()[:,:])

    my_avg = my_distrib.calc_averages_6d()
    my_rms = my_distrib.calc_rms_values_6d()

#    print('\nmy_avg =', my_avg[:])
#    print('\nmy_rms =', my_rms[:])

    for i in range(6):
        assert(my_avg[i] < 1.e-05)
        assert(stats6d.specify_significant_figures(my_rms[i], 3) == 1.0)

# test_unif_sphere()
