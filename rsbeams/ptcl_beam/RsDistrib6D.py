# -*- coding: utf-8 -*-
"""Generate a Gaussian or uniformly-filled 6D distribution.

Original code taken from RadTrack project, https://github.com/radiasoft/radtrack
:copyright: Copyright (c) 2013 RadiaBeam Technologies, LLC. All Rights Reserved.

Subsequent mods are due to RadiaSoft,
:copyright: Copyright (c) 2017 Radiasoft LLC. All Rights Reserved.

:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""
import math
import numpy
from rsbeams.ptcl_beam import RsTwiss2D
from rsbeams.ptcl_beam import RsPhaseSpace6D
from rsbeams.statistics import stats6d

class RsDistrib6D:
    """Generate a Gaussian or uniformly-filled 6D distribution."""

    def __init__(self, num_ptcls):
        # for testing purposes only
        if False:
            print ' '
            print ' ...in RsDistrib6D:__init__'
            print ' phase_space_6d object will be instantiated!'

        self.phase_space_6d = RsPhaseSpace6D.RsPhaseSpace6D(num_ptcls)
        self.phase_space_6d.check_array()

        # set some defaults
        self.maxRmsFactor = 5.0
        self.distrib_type = 'gaussian'
        return

    def get_phase_space_6d(self):
        return self.phase_space_6d

    def get_distrib_type(self):
        return self.distrib_type

    def set_distrib_type(self, distrib_type):
        if ( (distrib_type == 'uniform')  or
             (distrib_type == 'gaussian')):
            self.distrib_type = distrib_type
        else:
            message = '\n\nERROR --'
            message += '\n    distrib_type is specified as "' + distrib_type + '", which is not supported.'
            message += '\n    Only "uniform" and "gaussian" are allowed.'
            message += '\n'
            raise Exception(message)
        return

    def get_max_rms_fac(self):
        return self.max_rms_fac

    def set_max_rms_fac(self, max_rms_fac):
        # error handling of input data
        if (max_rms_fac > 0.0):
            self.max_rms_fac = max_rms_fac
        else:
            message = 'max_rms_fac = ' + str(max_rms_fac) + '; must be > 0.'
            raise Exception(message)
        return

    def make_unif_distrib(self):
        array_6d = self.phase_space_6d.get_array_6d()
        num_ptcls = self.phase_space_6d.get_num_ptcls()
        num_inside_circle = 0
        while (num_inside_circle < num_ptcls):
            test_x = 2. * numpy.random.uniform(0.0,1.0,1) - 1.
            test_y = 2. * numpy.random.uniform(0.0,1.0,1) - 1.
            test_z = 2. * numpy.random.uniform(0.0,1.0,1) - 1.
            test_sum = test_x**2 + test_y**2 + test_z**2

            if (test_sum < 1.):
                array_6d[0, num_inside_circle] = test_x
                array_6d[2, num_inside_circle] = test_y
                array_6d[4, num_inside_circle] = test_z
                num_inside_circle += 1

        num_inside_circle = 0
        while (num_inside_circle < num_ptcls):
            test_px = 2. * numpy.random.uniform(0.0,1.0,1) - 1.
            test_py = 2. * numpy.random.uniform(0.0,1.0,1) - 1.
            test_pz = 2. * numpy.random.uniform(0.0,1.0,1) - 1.
            test_sum = test_px**2 + test_py**2 + test_pz**2

            if (test_sum < 1.):
                array_6d[1, num_inside_circle] = test_px
                array_6d[3, num_inside_circle] = test_py
                array_6d[5, num_inside_circle] = test_pz
                num_inside_circle += 1

        return

    def make_gauss_distrib(self):
        array_6d = self.phase_space_6d.get_array_6d()
        num_ptcls = self.phase_space_6d.get_num_ptcls()
        for nLoop in range(6):
            num_inside_circle = 0
            while (num_inside_circle < num_ptcls):
                test_point = numpy.random.normal(0.0, 1.0, 1)

                if (test_point*test_point < self.max_rms_fac):
                    array_6d[nLoop, num_inside_circle] = test_point
                    num_inside_circle += 1
        return

    def init_phase_space_6d(self):
        if (self.distrib_type == 'uniform'):
            self.make_unif_distrib()
        elif (self.distrib_type == 'gaussian'):
            self.make_gauss_distrib()
        else:
            message = '\n\nERROR --'
            message += '\n    distrib_type is specified as "' + self.distrib_type + '", which is not supported.'
            message += '\n    Only "uniform" and "gaussian" are allowed.'
            message += '\n'
            raise Exception(message)
        return

    def clean_phase_space_6d(self):
        stats6d.sub_avg6d(self.phase_space_6d.get_array_6d())
        stats6d.rm_correlations6d(self.phase_space_6d.get_array_6d())
        stats6d.sub_avg6d(self.phase_space_6d.get_array_6d())
        stats6d.normalize_rms6d(self.phase_space_6d.get_array_6d())
        return

    def round_phase_space_6d(self):
        self.init_phase_space_6d()
        self.clean_phase_space_6d()
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
        for iLoop in range(3):
            ii = 2 * iLoop
            emitSQ = sigma[ii,ii]*sigma[ii+1,ii+1] - sigma[ii,ii+1]*sigma[ii+1,ii]

            if False:
                print ' '
                print ' num_ptcls = ', self.phase_space_6d.get_num_ptcls()
                q6 = self.phase_space_6d.get_array_6d()
                print ' 1st particle: ', q6[:,0]

            if False:
                print ' '
                print ' iLoop, ii = ', iLoop, ii
                print ' sigma[', ii,   ii,  '] = ', sigma[ii,  ii  ]
                print ' sigma[', ii+1, ii,  '] = ', sigma[ii+1,ii  ]
                print ' sigma[', ii,   ii+1,'] = ', sigma[ii,  ii+1]
                print ' sigma[', ii+1, ii+1,'] = ', sigma[ii+1,ii+1]

            if emitSQ <= 0.0:
                message  = 'Error -- \n\n'
                message += '  emitSQ = ' + str(emitSQ) + ' must be > zero!\n'
                message += '  ...in RsDistrib6D:calc_twiss6d()\n'
                message += '  iLoop, ii = ' + str(iLoop) + ', ' + str(ii) + '\n'
                raise Exception(message)

            emit_rms[iLoop]  =  math.sqrt(emitSQ)
            beta_rms[iLoop]  =  sigma[ii,ii]   / emit_rms[iLoop]
            alpha_rms[iLoop] = -sigma[ii,ii+1] / emit_rms[iLoop]

            if False:
                print ' '
                print ' alpha_rms, beta_rms, emit_rms = ', alpha_rms[iLoop], beta_rms[iLoop], emit_rms[iLoop]

        twiss6d['twissX'] = RsTwiss2D.RsTwiss2D(alpha_rms[0], beta_rms[0], emit_rms[0])
        twiss6d['twissY'] = RsTwiss2D.RsTwiss2D(alpha_rms[1], beta_rms[1], emit_rms[1])
        twiss6d['twissZ'] = RsTwiss2D.RsTwiss2D(alpha_rms[2], beta_rms[2], emit_rms[2])
        return

    def make_twiss_dist_6d(self,twiss6d, mean_p_ev):
        self.round_phase_space_6d()

        array_6d = self.phase_space_6d.get_array_6d()
        temp6D = array_6d.copy()

        ii = -1
        for iLoop in range(0,5,2):

            ii +=1
            if   ii==0: twissObject = twiss6d['twissX']
            elif ii==1: twissObject = twiss6d['twissY']
            elif ii==2: twissObject = twiss6d['twissZ']
            else:
                message = 'Error:  ii = ' + ii + ' -- not valid.'
                raise Exception(message)

            alphaII = twissObject.getAlphaRMS()
            betaII  = twissObject.getBetaRMS()
            gammaII = (1.0 + alphaII**2) / betaII

            if 0:
                print ' '
                print ' alpha, beta, gamma[', ii, '] = ', alphaII, betaII, gammaII

            gMinusB = gammaII - betaII
            rootFac = math.sqrt(gMinusB**2 + 4.0*alphaII**2)

            if 0:
                print ' gMinusB, rootFac[', ii, '] = ', gMinusB, rootFac

            if gMinusB >= 0.0:
                fac  = math.sqrt(0.5*(gammaII+betaII-rootFac))
                fInv = math.sqrt(0.5*(gammaII+betaII+rootFac))
            else:
                fac  = math.sqrt(0.5*(gammaII+betaII+rootFac))
                fInv = math.sqrt(0.5*(gammaII+betaII-rootFac))

            if 0:
                print ' fac, fInv[', ii, '] = ', fac, fInv

            if alphaII == 0.0:
                sinPhi = 0.0
                cosPhi = 1.0
            else:
                sinPhi = math.sqrt(0.5*(1.-math.fabs(gMinusB)/rootFac))
                cosPhi = math.sqrt(0.5*(1.+math.fabs(gMinusB)/rootFac))

            if alphaII*gMinusB < 0.0: sinPhi = -sinPhi

            rootFac = math.sqrt(twissObject.getEmitRMS())

            if 0:
                print ' sinPhi, cosPhi, rootFac[', ii, '] = ', sinPhi, cosPhi, rootFac

            for nLoop in range(self.phase_space_6d.getNumParticles()):
                array_6d[iLoop  ,nLoop] = rootFac*(fac *cosPhi*temp6D[iLoop,  nLoop] - \
                                                  fInv*sinPhi*temp6D[iLoop+1,nLoop])
                array_6d[iLoop+1,nLoop] = rootFac*(fac *sinPhi*temp6D[iLoop,  nLoop] + \
                                                  fInv*cosPhi*temp6D[iLoop+1,nLoop])
        self.multiply_distrib_component(mean_p_ev, 5)
        self.offset_distrib_component(mean_p_ev, 5)

    def offset_distrib_component(self,offset,index):
        if index < 0 or index > 5:
            message = 'ERROR!  index is out of range: ' + str(index)
            raise Exception(message)

        array_6d = self.phase_space_6d.get_array_6d()
        array_6d[index,:] += offset

        return

    def multiply_distrib_component(self,factor,index):
        if index < 0 or index > 5:
            message = 'ERROR!  index is out of range: ' + str(index)
            raise Exception(message)

        array_6d = self.phase_space_6d.get_array_6d()
        array_6d[index,:] *= factor

        return
