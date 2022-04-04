# -*- coding: utf-8 -*-
"""DEPRECATED use rsmath.const or scipy

:copyright: Copyright (c) 2017 Radiasoft LLC. All Rights Reserved.
:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""
import math
import scipy.constants as const

# DEPRECATED use rsmath.rsconst or scipy
C_to_STATC = const.c * 10.
EV_to_ERG = 1.602e-12
KG_to_EV = c_SQ / e
MKS_factor = 1./(4.*math.pi*const.epsilon_0)
RT_2 = math.sqrt(2.)
RT_2_OVER_PI = math.sqrt(2/math.pi)
RT_TWO_PI = math.sqrt(2*math.pi)
TWO_PI = 2 * math.pi
c = const.c
c_CGS = const.c * 100
c_INV  = 1./const.c
c_SQ = const.c**2
e = const.e
e_CGS = const.e * C_to_STATC
epsilon_0 = const.epsilon_0
m_e = const.m_e
m_e_CGS = const.m_e * 1000.
m_e_EV = const.m_e * KG_to_EV
m_p = const.m_p
m_p_CGS = const.m_p * 1000.
m_p_EV = const.m_p * KG_to_EV
mu_0 = const.mu_0
pi = math.pi