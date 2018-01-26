from __future__ import absolute_import, division, print_function, unicode_literals
import pytest
import numpy
from rsbeams.rsstats import stats6d

def test_stats_01():
    array_6d = numpy.zeros(36).reshape(6,6)
    array_6d[0,0] = 1.
    array_6d[1,0] = 2.
    array_6d[2,0] = 3.
    array_6d[3,0] = 4.
    array_6d[4,0] = 5.
    array_6d[5,0] = 6.
#    print('array_6d = ')
#    print(array_6d[:,:])

    avg_test = numpy.zeros(6)
    for i in range(6):
        avg_test[i] = (i+1)/6.
#    print('avg_test = ')
#    print(avg_test[:])

    avg_6d = stats6d.calc_avg6d(array_6d)
#    print('avg_6d = ')
#    print(avg_6d[:])

#    for i in range(6):
#        print('avg_6d[',i,']   = ',avg_6d[i])
#        print('avg_test[',i,'] = ',avg_test[i])

    for i in range(6):
        assert stats6d.specify_significant_figures(avg_6d[i], 4) == \
            stats6d.specify_significant_figures(avg_test[i], 4)

#test_stats_01()
