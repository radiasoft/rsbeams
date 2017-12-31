from __future__ import absolute_import, division, print_function, unicode_literals
import pytest
import numpy
from rsbeams.statistics import stats6d
from rsbeams.ptcl_beam import RsPhaseSpace6D

def test_phase_space_01():
    # initialization creates a 6 x N array of zeros
    my_n = 6
    my_phase_space = RsPhaseSpace6D.RsPhaseSpace6D(my_n)
    test_name = my_phase_space.get_file_name()
    assert(test_name == 'ps6d')
    test_array = my_phase_space.get_array_6d()
    my_shape = numpy.shape(test_array)
    assert(len(my_shape) == 2)
    assert(my_shape[0] == 6)
    assert(my_shape[1] == my_n)

    my_array = numpy.zeros(36).reshape(6,6)
    my_array[0,0] = 1.
    my_array[1,0] = 2.
    my_array[2,0] = 3.
    my_array[3,0] = 4.
    my_array[4,0] = 5.
    my_array[5,0] = 6.
    my_phase_space.set_array_6d(my_array)
    my_ps_0 = my_phase_space.get_array_x()
    assert(my_ps_0[0] == my_array[0,0])
    my_ps_1 = my_phase_space.get_array_xp()
    assert(my_ps_1[0] == my_array[1,0])
    my_ps_2 = my_phase_space.get_array_y()
    assert(my_ps_2[0] == my_array[2,0])
    my_ps_3 = my_phase_space.get_array_yp()
    assert(my_ps_3[0] == my_array[3,0])
    my_ps_4 = my_phase_space.get_array_s()
    assert(my_ps_4[0] == my_array[4,0])
    my_ps_5 = my_phase_space.get_array_dp()
    assert(my_ps_5[0] == my_array[5,0])

# test_phase_space_01()
