# -*- coding: utf-8 -*-
"""Encapsulation of statistics for a 6D particle distribution.

Original code taken from RadTrack project, https://github.com/radiasoft/radtrack
:copyright: Copyright (c) 2013 RadiaBeam Technologies, LLC. All Rights Reserved.

Subsequent mods are due to RadiaSoft,
:copyright: Copyright (c) 2017 Radiasoft LLC. All Rights Reserved.

:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""
import math
import numpy
import scipy

def calc_avg6d(array6D):
    return scipy.average(array6D, axis=1)

def sub_avg6d(array6D):
    avg6d = calc_avg6d(array6D)
    for nLoop in range(6):
        array6D[nLoop,:] -= avg6d[nLoop]

def calc_variance6d(array6D):
    npoints = array6D.shape[1]
    avgIter = array6D[:,0]
    varIter = scipy.zeros(6)
    for nLoop in range(1,npoints):
        tmpData  = array6D[:,nLoop]
        tmpIter  = (tmpData - avgIter)
        avgIter += (tmpData - avgIter) / (nLoop+1)
        varIter += (tmpData - avgIter) * tmpIter
    return varIter / npoints

def calc_rms6d(array6D):
    return scipy.sqrt(calc_variance6d(array6D))

def normalize_rms6d(array6D):
    invRmsValues6D = 1. / calc_rms6d(array6D)
    for nLoop in range(6):
        array6D[nLoop,:] *= invRmsValues6D[nLoop]

def calc_min6d(array6D):
    return numpy.min(array6D, axis=1)

def calc_max6d(array6D):
    return numpy.max(array6D, axis=1)

def calc_correlations6d(array6D):
    npoints = array6D.shape[1]
    avg6d = calc_avg6d(array6D)
    variance6D = calc_variance6d(array6D)
    correlations6D = scipy.zeros(6*6).reshape(6,6)
    for iLoop in range(6):
        for jLoop in range(6):
            if iLoop == jLoop:
                correlations6D[iLoop, jLoop] = variance6D[iLoop]
            else:
                for nLoop in range(npoints):
                    correlations6D[iLoop, jLoop] += \
                        (array6D[iLoop,nLoop] - avg6d[iLoop]) * \
                        (array6D[jLoop,nLoop] - avg6d[jLoop])
                correlations6D[iLoop, jLoop] /= npoints
    return correlations6D

def rm_correlations6d(array6D):
    npoints = array6D.shape[1]
    sigmaM = calc_correlations6d(array6D)
    eigVals, eigVecs = jacobi_eigen_solver6d(sigmaM)

    verboseCheck = 0
    if verboseCheck == 1:
        print 'eigVals = ', eigVals

    temp6D = array6D.copy()
    for iLoop in range(6):
        for nLoop in range(npoints): array6D[iLoop,nLoop] = 0.0

    for iLoop in range(6):
        for jLoop in range(6):
            for nLoop in range(npoints):
                array6D[iLoop,nLoop] += eigVecs[jLoop,iLoop] * temp6D[jLoop,nLoop]

def jacobi_eigen_solver6d(sigma6D):
    # Setup
    eVecs=scipy.zeros(36).reshape(6,6)
    for ip in range(6): eVecs[ip,ip]=1.0

    bTemp=scipy.zeros(6)
    zTemp=scipy.zeros(6)
    eVals=scipy.zeros(6)
    for ip in range(6):
        bTemp[ip]=sigma6D[ip,ip]
        eVals[ip]=sigma6D[ip,ip]

    # Top of the master loop
    numRotations = 0
    for nMaster in range(50):
        sm = 0.0
        for ip in range(5):
            for iq in range(ip+1,6):
                sm += math.fabs(sigma6D[ip,iq])

        # Check for convergence
        if sm == 0.0:
            return eVals, eVecs   # Success!

        # Convergence failed, so reset threshold
        if nMaster<3:
            threshold=0.2*sm/36.
        else:
            threshold=0.0

        # Next iteration
        for ip in range(5):
            for iq in range(ip+1,6):
                gScal=100.*math.fabs(sigma6D[ip,iq])
                if nMaster>3 and math.fabs(float(eVals[ip])+gScal)==math.fabs(eVals[ip]) \
                             and math.fabs(float(eVals[iq])+gScal)==math.fabs(eVals[iq]):
                    sigma6D[ip,iq]=0.0
                elif math.fabs(sigma6D[ip,iq])>threshold:
                    hScal=float(eVals[iq])-float(eVals[ip])
                    if math.fabs(hScal)+gScal==math.fabs(hScal):
                        tScal=float(sigma6D[ip,iq])/hScal
                    else:
                        theta=0.5*hScal/float(sigma6D[ip,iq])
                        tScal=1.0/(math.fabs(theta)+math.sqrt(1.0+theta**2))
                        if theta<0.: tScal*=-1.0
                    cTemp=1.0/math.sqrt(1.0+tScal**2)
                    sTemp=tScal*cTemp
                    tau=sTemp/(1.0+cTemp)
                    hScal=tScal*float(sigma6D[ip,iq])
                    zTemp[ip]-=hScal
                    zTemp[iq]+=hScal
                    eVals[ip]-=hScal
                    eVals[iq]+=hScal
                    sigma6D[ip,iq]=0.0
                    for jLoop in range(ip):
                        gScal=sigma6D[jLoop,ip]
                        hScal=sigma6D[jLoop,iq]
                        sigma6D[jLoop,ip]=gScal-sTemp*(hScal+gScal*tau)
                        sigma6D[jLoop,iq]=hScal+sTemp*(gScal-hScal*tau)
                    for jLoop in range(ip+1,iq):
                        gScal=sigma6D[ip,jLoop]
                        hScal=sigma6D[jLoop,iq]
                        sigma6D[ip,jLoop]=gScal-sTemp*(hScal+gScal*tau)
                        sigma6D[jLoop,iq]=hScal+sTemp*(gScal-hScal*tau)
                    for jLoop in range(iq+1,6):
                        gScal=sigma6D[ip,jLoop]
                        hScal=sigma6D[iq,jLoop]
                        sigma6D[ip,jLoop]=gScal-sTemp*(hScal+gScal*tau)
                        sigma6D[iq,jLoop]=hScal+sTemp*(gScal-hScal*tau)
                    for jLoop in range(6):
                        gScal=eVecs[jLoop,ip]
                        hScal=eVecs[jLoop,iq]
                        eVecs[jLoop,ip]=gScal-sTemp*(hScal+gScal*tau)
                        eVecs[jLoop,iq]=hScal+sTemp*(gScal-hScal*tau)
                    numRotations+=1

        # Collect results before checking again for convergence
        for ip in range(6):
            bTemp[ip]+=zTemp[ip]
            eVals[ip] =bTemp[ip]
            zTemp[ip] =0.0

    # No convergence after 50 iterations, so give up!
    raise Exception("Too many iterations in routine jacobi_eigen_solver6d")

def specify_significant_figures(float_value, num_sig_figs):
    """Round float_value to num_sig_figs significant figures."""
    try:
        # Find my_mult, my_exp such that float_value = my_mult*10^my_exp
        #   such that  1. <= my_mult < 10.
        my_exp = math.floor(math.log10(abs(float_value)))
        my_mult = float_value/(10**my_exp)
        return round(my_mult, num_sig_figs-1)*(10**my_exp)
    except ValueError:
        return 0
