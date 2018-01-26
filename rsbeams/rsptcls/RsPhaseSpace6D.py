# -*- coding: utf-8 -*-
"""Data holding class for 6D phase space distributions.

Original code taken from RadTrack project, https://github.com/radiasoft/radtrack
:copyright: Copyright (c) 2013 RadiaBeam Technologies, LLC. All Rights Reserved.

Subsequent mods are due to RadiaSoft,
:copyright: Copyright (c) 2017 Radiasoft LLC. All Rights Reserved.

:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""

import numpy

class RsPhaseSpace6D:
    """Data holding class for 6D phase space distributions."""

    def __init__(self, num_ptcls):
        # set the number of macro-particles in the bunch
        # initialize the 6D phase space array
        self.set_num_ptcls(num_ptcls)

        # specify default file name
        self.file_name = 'ps6d'

        return

    def get_num_ptcls(self):
        return self.num_ptcls

    def set_num_ptcls(self, num_ptcls):
        if (num_ptcls > 0):
            self.num_ptcls = num_ptcls
            self.initialize_array()
        else:
            message = 'ERROR - num_ptcls <= 0: ' + str(num_ptcls)
            raise Exception(message)
        return

    def initialize_array(self):
        self.array_6d = numpy.zeros(6*self.num_ptcls).reshape(6, self.num_ptcls)
        return

    def get_file_name(self):
        return self.file_name

    def set_file_name(self, file_name):
        # error handling of input data
        if (len(file_name) > 0):
            self.file_name = file_name
        else:
            message = 'file_name must have length > 0'
            raise Exception(message)
        return

    def write_array(self):
        numpy.savez(self.file_name+'.npz', a=self.array_6d)
        return

    def read_array(self):
        dataObject = numpy.load(self.file_name + '.npz')
        self.array_6d = dataObject['a']
        self.check_array()
        return

    def check_array(self):
        array_shape = numpy.shape(self.array_6d)

        if len(array_shape) != 2:
            message = 'Dimensionality of array_6d matrix = ' + str(len(array_shape)) + ' -- must be 2.'
            raise Exception(message)
        else:
            if array_shape[0] != 6:
                message = 'Dimensionality of array_6d = ' + str(array_shape[0]) + ' -- must be 6.'
                raise Exception(message)
            if (array_shape[1] < 1):
                message = 'Number of array_6d points = ' + str(array_shape[1]) + ' -- must be > 0'
                raise Exception(message)
            else:
                self.num_ptcls = array_shape[1]

        return

    def get_array_6d(self):
        return self.array_6d

    def set_array_6d(self,array_6d):
        self.array_6d = array_6d
        return

    def get_array_x(self):
        return self.array_6d[0,:]

    def get_array_xp(self):
        return self.array_6d[1,:]

    def get_array_y(self):
        return self.array_6d[2,:]

    def get_array_yp(self):
        return self.array_6d[3,:]

    def get_array_s(self):
        return self.array_6d[4,:]

    def get_array_dp(self):
        return self.array_6d[5,:]
