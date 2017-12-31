# -*- coding: utf-8 -*-
"""Wrapper for scipy.constants so we can add more constants

We obey the PEP-8 convention that constants should be all caps with underscores,
except when we are following scipy or modifying a scipy constant.

:copyright: Copyright (c) 2017 Radiasoft LLC. All Rights Reserved.
:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""
import math
import scipy.constants

# obviously ok to use math.pi everywhere; here for completeness
pi = math.pi
RT_2 = math.sqrt(2.)
TWO_PI = 2 * math.pi
RT_TWO_PI = math.sqrt(2*math.pi)
RT_2_OVER_PI = math.sqrt(2/math.pi)

c = const.c                  # speed of light [m/s]
c_SQ = const.c**2
c_INV  = 1./const.c

mu_0 = const.mu_0            # permeability of free space
epsilon_0 = const.epsilon_0  # permittivity of free space

m_p = const.m_p              # proton mass [kg]
m_e = const.m_e              # electron mass [kg]
e = const.e                  # fundamental electric charge [C] (positive)

KG_to_EV = c_SQ / e          # convert mass [kg] to effective energy [eV]
EV_to_ERG = 1.602e-12        # convert energy from eV to erg
C_to_STATC = const.c * 10.   # convert charge from coulombs to statcoulombs

e_CGS = const.e * C_to_STATC
c_CGS = const.c * 100
m_e_CGS = const.m_e * 1000.
m_p_CGS = const.m_p * 1000.
m_e_EV = const.m_e * KG_to_EV
m_p_EV = const.m_p * KG_to_EV
