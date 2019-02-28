# -*- coding: utf-8 -*-
"""Calculations from 2017 PRAB article by Lebedev, Burov and Nagaitsev (LBN)

:copyright: Copyright (c) 2019 Radiasoft LLC. All Rights Reserved.
:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""
# import the usual suspects
from __future__ import print_function
import math
import numpy as np
from rsbeams.rsphysics import rsconst


# Specify default values of some physical quantities...
#   These are motivated by upcoming PWFA experiments at FACET II

# number density of pre-ionized electron plasma
n_pe_cgs = 4.e16             # [cm^-3]
n_pe     = n_pe_cgs * 1.e6   #  [m^-3]

# electron plasma frequency, wavenumber, wavelength
om_pe     = np.sqrt(n_pe*rsconst.e**2
            / (rsconst.m_e*rsconst.epsilon_0) )  # [rad/s]
k_pe      = om_pe/constants.c                    # [rad/m]
lambda_pe = 2.*np.pi/k_pe                        # [m]
lambda_pe_microns = lambda_pe*1.e6               # [microns]

# -------------------
# The PWFA drive beam, with an assumed Gaussian distribution
# -------------------

# RMS radius
# beam_rms_r = 3.65e-6   # [m]   FACET II params...?
beam_rms_r = 0.4 / k_pe  # [m]   satisfies resonance condition

# RMS length
# beam_rms_z = 12.77e-6  # [m]   FACET II params...?
beam_rms_z = 1.2 / k_pe  # [m]   satisfies resonance condition

# total charge in the electron drive beam
# beam_tot_q = 1.e10*math.abs(rsconst.e)    # [C]   FACET II params...?
beam_tot_q = 3.e-9                          # [C]   3 nC

# peak number density of the drive beam
beam_rms_v = math.pi*beam_rms_z*beam_rms_r**2
beam_rms_n = beam_tot_q/math.abs(rsconst.e)/beam_rms_v

# large density ratio is required for blowout regime
dens_ratio = beam_rms_n/n_pe    # should be >> unity

# relativistic gamma factor
beam_gamma = 1.957e+4    # FACET II params...?

# -------------------
# The accelerated "witness" beam, also assumed Gaussian
# -------------------

wb_rms_r = 0.012 * lambda_pe  # [m] RMS radius
wb_rms_z = 0.036 * lambda_pe  # [m] RMS length
wb_tot_q = 0.100e-9           # [C] total charge of 100 pC
wb_gamma = 100                # relativistic gamma factor
wb_trail = 0.900 * lambda_pe  # [m] trailing distance behind center of drive beam

def calc_Ez_on_axis(n_pe, rb, drb_dxi):
    """
    Calculate the longitudinal electric field in a plasma bubble.
    Valid in "strong bubble regime", where rb_max*k_pe >> 1.
    The calculation is local; it depends on bubble radius and slope.
    The slope of rb is calculated wrt xi, where xi=ct-z

    Args:
        n_pe:    number density of the electron plasma
        rb:      the local bubble radius
        drb_dxi: longitudinal derivative of the bubble radius

    Returns:
        Ez: longitudinal electric field along the axis
    """
    Ez = -2.*math.pi*n_pe*math.abs(rsconst.e)*rb*drb_dxi
    return Ez

def calc_drb_dxi_no_beam(rb, rb_max):
    """
    Calculate the longitudinal slope of the local bubble radius.
    Valid in "strong bubble regime", where rb_max*k_pe >> 1.
    Valid only in back of bubble, where there is no drive beam.
    The slope of rb is calculated wrt xi, where xi=ct-z

    Args:
        rb:     the local bubble radius
        rb_max: maximum value of the bubble radius

    Returns:
        drb_dxi: longitudinal derivative of the bubble radius
    """
    # there is ambiguity in the sign, which needs to be resolved
    drb_dxi = math.sqrt((pow(rb_max/r_b,4)-1.)/2.)
    return drb_dxi

def calc_Ez_on_axis_no_beam(n_pe, rb, drb_dxi):
    """
    Calculate the longitudinal electric field in a plasma bubble.
    Valid in "strong bubble regime", where rb_max*k_pe >> 1.
    Valid only in back of bubble, where there is no drive beam.
    The calculation is local; it depends on bubble radius and slope.

    Args:
        n_pe:    number density of the electron plasma
        rb:      the local bubble radius
        drb_dxi: longitudinal derivative of the bubble radius

    Returns:
        Ez: longitudinal electric field along the axis
    """
    Ez = -2.*math.pi*n_pe*math.abs(rsconst.e)*rb*drb_dxi
    return Ez
