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

class RsPlotPhaseSpace6D:
    """Plot phase space points from RsPhaseSpace6D class"""

    def __init__(self, phaseSpace6D):

        # error checking
        phaseSpace6D.check_array()
        tmpPoints = phaseSpace6D.get_num_ptcls()
        if (tmpPoints > 1):
            self.npoints = tmpPoints
        else:
            message = 'ERROR - phaseSpace6D.getNumParticles() < 2: ' + str(tmpPoints)
            raise Exception(message)

        self.data6D = phaseSpace6D.get_array_6d()

        self.label=numpy.array(['x [m]', 'px', 'y [m]', 'py', 'z [m]', 'pz'])
        self.title='no plot title specified'
        self.figNum = 0
        return

    def setTitle(self,title):
        self.title = title
        return

    def getTitle(self):
        return self.title

    def setLabel(self,index,label):
        self.label[index] = label
        return

    def getLabel(self,index):
        return self.label[index]

    def getData6D(self):
        return self.data6D

    def setData6D(self,data6D):
        self.data6D = data6D
        return

    def plotData6D(self,hIndex,vIndex):
        hArray = self.data6D[hIndex,:]
        vArray = self.data6D[vIndex,:]

        hMin = min(hArray)
        hMax = max(hArray)
        if -hMin > hMax:
            hMax = math.fabs(hMin)
        else:
            hMin = -hMax

        vMin = min(vArray)
        vMax = max(vArray)
        if -vMin > vMax:
            vMax = math.fabs(vMin)
        else:
            vMin = -vMax

        self.figNum += 1
        pyplot.figure(self.figNum)
        pyplot.scatter(hArray, vArray, marker=',',s=1, c='k')
        pyplot.axis([hMin, hMax, vMin, vMax])

        pyplot.xlabel(self.label[hIndex])
        pyplot.ylabel(self.label[vIndex])
        pyplot.title(self.title)

        return

    def showPlots(self):
        pyplot.show()
        return
