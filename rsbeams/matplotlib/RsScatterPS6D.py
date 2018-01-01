# -*- coding: utf-8 -*-
"""Plot phase space points from RsPhaseSpace6D class

Original code taken from RadTrack project, https://github.com/radiasoft/radtrack
:copyright: Copyright (c) 2014 RadiaBeam Technologies, LLC. All Rights Reserved.

Subsequent mods are due to RadiaSoft,
:copyright: Copyright (c) 2017 Radiasoft LLC. All Rights Reserved.

:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""
import math
import numpy
from matplotlib import pyplot

class RsScatterPS6D:
    """Enable scatter plots of points from RsPhaseSpace6D class"""

    def __init__(self, phaseSpace6D):

        # error checking
        phaseSpace6D.check_array()
        tmpPoints = phaseSpace6D.get_num_ptcls()
        if (tmpPoints > 1):
            self.npoints = tmpPoints
        else:
            message = 'ERROR - phaseSpace6D.getNumParticles() < 2: ' + str(tmpPoints)
            raise Exception(message)

        self.data6d = phaseSpace6D.get_array_6d()

        self.label=numpy.array(['x [m]', 'px/p0', 'y [m]', 'py/p0', 's [m]', 'dp/p0'])
        self.title='no plot title specified'
        self.figNum = 0
        return

    def set_title(self,title):
        self.title = title
        return

    def get_title(self):
        return self.title

    def set_label(self,index,label):
        self.label[index] = label
        return

    def get_label(self,index):
        return self.label[index]

    def get_data6d(self):
        return self.data6d

    def set_data6d(self,data6d):
        self.data6d = data6d
        return

    def plot_data6d(self,h_index,v_index):
        h_array = self.data6d[h_index,:]
        v_array = self.data6d[v_index,:]

        h_min = min(h_array)
        h_max = max(h_array)
        if -h_min > h_max:
            h_max = math.fabs(h_min)
        else:
            h_min = -h_max

        v_min = min(v_array)
        v_max = max(v_array)
        if -v_min > v_max:
            v_max = math.fabs(v_min)
        else:
            v_min = -v_max

        self.figNum += 1
        pyplot.figure(self.figNum)
        pyplot.scatter(h_array, v_array, marker=',',s=1, c='k')
        pyplot.axis([h_min, h_max, v_min, v_max])

        pyplot.xlabel(self.label[h_index])
        pyplot.ylabel(self.label[v_index])
        pyplot.title(self.title)

        return

    def show_plots(self):
        pyplot.show()
        return

    def clear_plots(self):
        pyplot.clf()
        pyplot.cla()
        pyplot.close()
        return

    def save_plots(self,file_name):
        pyplot.savefig(file_name, bbox_inches='tight')
        pyplot.clf()
        pyplot.cla()
        return
