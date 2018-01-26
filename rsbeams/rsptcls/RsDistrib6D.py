# -*- coding: utf-8 -*-
"""Generate a Gaussian or uniformly-filled 6D distribution.

Original code taken from RadTrack project, https://github.com/radiasoft/radtrack
:copyright: Copyright (c) 2013 RadiaBeam Technologies, LLC. All Rights Reserved.

Subsequent mods are due to RadiaSoft,
:copyright: Copyright (c) 2017 Radiasoft LLC. All Rights Reserved.

:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import math
import numpy
from rsbeams.rsptcls import RsTwiss2D
from rsbeams.rsptcls import RsPhaseSpace6D
from rsbeams.rsstats import stats6d

class RsDistrib6D:
    """Generate a Gaussian or uniformly-filled 6D distribution."""

    def __init__(self, num_ptcls, distrib_type, max_rms_fac):

        self.phase_space_6d = RsPhaseSpace6D.RsPhaseSpace6D(num_ptcls)
        self.phase_space_6d.check_array()

        if ( (distrib_type != 'uniform') and
             (distrib_type != 'gaussian') ):
            message = '\n\nERROR --'
            message += '\n    distrib_type is specified as "' + self.distrib_type + '", which is not supported.'
            message += '\n    Only "uniform" and "gaussian" are allowed.'
            message += '\n'
            raise Exception(message)
        else:
            self.distrib_type = distrib_type

        if (max_rms_fac > 0.0):
            self.max_rms_fac = max_rms_fac
        else:
            message = 'max_rms_fac = ' + str(max_rms_fac) + '; must be > 0.'
            raise Exception(message)

        if (self.distrib_type == 'uniform'):
            self.make_unif_distrib()

        if (self.distrib_type == 'gaussian'):
            self.make_gauss_distrib()

        self.clean_phase_space_6d()
        return

    def get_phase_space_6d(self):
        return self.phase_space_6d

    def get_distrib_type(self):
        return self.distrib_type

    def get_max_rms_fac(self):
        return self.max_rms_fac

    def make_unif_distrib(self):
        array6d = self.phase_space_6d.get_array_6d()
        num_ptcls = self.phase_space_6d.get_num_ptcls()
        num_inside_circle = 0
        while (num_inside_circle < num_ptcls):
            test_x = 2. * numpy.random.uniform(0.0,1.0,1) - 1.
            test_y = 2. * numpy.random.uniform(0.0,1.0,1) - 1.
            test_z = 2. * numpy.random.uniform(0.0,1.0,1) - 1.
            test_sum = test_x**2 + test_y**2 + test_z**2

            if (test_sum < 1.):
                array6d[0, num_inside_circle] = test_x
                array6d[2, num_inside_circle] = test_y
                array6d[4, num_inside_circle] = test_z
                num_inside_circle += 1

        num_inside_circle = 0
        while (num_inside_circle < num_ptcls):
            test_px = 2. * numpy.random.uniform(0.0,1.0,1) - 1.
            test_py = 2. * numpy.random.uniform(0.0,1.0,1) - 1.
            test_pz = 2. * numpy.random.uniform(0.0,1.0,1) - 1.
            test_sum = test_px**2 + test_py**2 + test_pz**2

            if (test_sum < 1.):
                array6d[1, num_inside_circle] = test_px
                array6d[3, num_inside_circle] = test_py
                array6d[5, num_inside_circle] = test_pz
                num_inside_circle += 1

        return

    def make_gauss_distrib(self):
        array6d = self.phase_space_6d.get_array_6d()
        num_ptcls = self.phase_space_6d.get_num_ptcls()
        for nLoop in range(6):
            num_inside_circle = 0
            while (num_inside_circle < num_ptcls):
                test_point = numpy.random.normal(0.0, 1.0, 1)

                if (test_point*test_point < self.max_rms_fac):
                    array6d[nLoop, num_inside_circle] = test_point
                    num_inside_circle += 1
        return

    def clean_phase_space_6d(self):
        stats6d.sub_avg6d(self.phase_space_6d.get_array_6d())
        stats6d.rm_correlations6d(self.phase_space_6d.get_array_6d())
        stats6d.sub_avg6d(self.phase_space_6d.get_array_6d())
        stats6d.normalize_rms6d(self.phase_space_6d.get_array_6d())
        return

    def calc_averages_6d(self):
        averages = stats6d.calc_avg6d(self.phase_space_6d.get_array_6d())
        return averages

    def calc_rms_values_6d(self):
        rmsValues = stats6d.calc_rms6d(self.phase_space_6d.get_array_6d())
        return rmsValues

    def calc_twiss6d(self,twiss6d):
        alpha_rms = numpy.zeros(3)
        beta_rms  = numpy.zeros(3)
        emit_rms  = numpy.zeros(3)

        sigma = stats6d.calc_correlations6d(self.phase_space_6d.get_array_6d())
        for i_loop in range(3):
            ii = 2 * i_loop
            emitSQ = sigma[ii,ii]*sigma[ii+1,ii+1] - sigma[ii,ii+1]*sigma[ii+1,ii]

            if False:
                print('\n num_ptcls = ', self.phase_space_6d.get_num_ptcls())
                q6 = self.phase_space_6d.get_array_6d()
                print(' 1st particle: ', q6[:,0])

            if False:
                print('\n i_loop, ii = ', i_loop, ii)
                print(' sigma[', ii,   ii,  '] = ', sigma[ii,  ii  ])
                print(' sigma[', ii+1, ii,  '] = ', sigma[ii+1,ii  ])
                print(' sigma[', ii,   ii+1,'] = ', sigma[ii,  ii+1])
                print(' sigma[', ii+1, ii+1,'] = ', sigma[ii+1,ii+1])

            if emitSQ <= 0.0:
                message  = 'Error -- \n\n'
                message += '  emitSQ = ' + str(emitSQ) + ' must be > zero!\n'
                message += '  ...in RsDistrib6D:calc_twiss6d()\n'
                message += '  i_loop, ii = ' + str(i_loop) + ', ' + str(ii) + '\n'
                raise Exception(message)

            emit_rms[i_loop]  =  math.sqrt(emitSQ)
            beta_rms[i_loop]  =  sigma[ii,ii]   / emit_rms[i_loop]
            alpha_rms[i_loop] = -sigma[ii,ii+1] / emit_rms[i_loop]

            if False:
                print('\n alpha_rms, beta_rms, emit_rms = ', \
                    alpha_rms[i_loop], beta_rms[i_loop], emit_rms[i_loop])

        twiss6d['twiss_x'] = RsTwiss2D.RsTwiss2D(alpha_rms[0], beta_rms[0], emit_rms[0])
        twiss6d['twiss_y'] = RsTwiss2D.RsTwiss2D(alpha_rms[1], beta_rms[1], emit_rms[1])
        twiss6d['twiss_z'] = RsTwiss2D.RsTwiss2D(alpha_rms[2], beta_rms[2], emit_rms[2])
        return

    def make_twiss_dist_6d(self,twiss6d, mean_p_ev):

        array6d = self.phase_space_6d.get_array_6d()
        temp6D = array6d.copy()

        ii = -1
        for i_loop in range(0,5,2):

            ii +=1
            if   ii==0: twissObject = twiss6d['twiss_x']
            elif ii==1: twissObject = twiss6d['twiss_y']
            elif ii==2: twissObject = twiss6d['twiss_z']
            else:
                message = 'Error:  ii = ' + ii + ' -- not valid.'
                raise Exception(message)

            alphaII = twissObject.get_alpha_rms()
            betaII  = twissObject.get_beta_rms()
            gammaII = (1.0 + alphaII**2) / betaII

            if 0:
                print('\n alpha, beta, gamma[', ii, '] = ', alphaII, betaII, gammaII)

            gMinusB = gammaII - betaII
            rt_fac = math.sqrt(gMinusB**2 + 4.0*alphaII**2)

            if 0:
                print(' gMinusB, rt_fac[', ii, '] = ', gMinusB, rt_fac)

            if gMinusB >= 0.0:
                fac  = math.sqrt(0.5*(gammaII+betaII-rt_fac))
                f_inv = math.sqrt(0.5*(gammaII+betaII+rt_fac))
            else:
                fac  = math.sqrt(0.5*(gammaII+betaII+rt_fac))
                f_inv = math.sqrt(0.5*(gammaII+betaII-rt_fac))

            if 0:
                print(' fac, f_inv[', ii, '] = ', fac, f_inv)

            if alphaII == 0.0:
                sin_phi = 0.0
                cos_phi = 1.0
            else:
                sin_phi = math.sqrt(0.5*(1.-math.fabs(gMinusB)/rt_fac))
                cos_phi = math.sqrt(0.5*(1.+math.fabs(gMinusB)/rt_fac))

            if alphaII*gMinusB < 0.0: sin_phi = -sin_phi

            rt_fac = math.sqrt(twissObject.get_emit_rms())

            if 0:
                print(' sin_phi, cos_phi, rt_fac[', ii, '] = ', sin_phi, cos_phi, rt_fac)

            for nLoop in range(self.phase_space_6d.get_num_ptcls()):
                array6d[i_loop  ,nLoop] = rt_fac*(fac * cos_phi*temp6D[i_loop,  nLoop] - \
                                                  f_inv*sin_phi*temp6D[i_loop+1,nLoop])
                array6d[i_loop+1,nLoop] = rt_fac*(fac * sin_phi*temp6D[i_loop,  nLoop] + \
                                                  f_inv*cos_phi*temp6D[i_loop+1,nLoop])
#        self.multiply_component(mean_p_ev, 5)
#        self.offset_component(mean_p_ev, 5)

        return

    def offset_component(self,offset,index):
        if index < 0 or index > 5:
            message = 'ERROR!  index is out of range: ' + str(index)
            raise Exception(message)

        array6d = self.phase_space_6d.get_array_6d()
        array6d[index,:] += offset

        return

    def multiply_component(self,factor,index):
        if index < 0 or index > 5:
            message = 'ERROR!  index is out of range: ' + str(index)
            raise Exception(message)

        array6d = self.phase_space_6d.get_array_6d()
        array6d[index,:] *= factor

        return
