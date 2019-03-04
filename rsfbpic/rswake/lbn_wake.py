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

def calc_rb_max(n_pe, beam_tot_z, beam_num_ptcl):
    """
    Calculate the maximum radius of the plasma bubble

    Valid in "strong bubble regime", where rb_max*k_pe >> 1.
    Args:
        n_pe:           number density of the electron plasma
        beam_tot_z:     total length of the drive beam
        beam_num_ptcl:  number of e- in the drive beam
    Returns:
        rb_max: maximum radius of the plasma bubble
    """
    # from Eq. (12) of LBN2017
    rb_max = pow(2,7/8.)*pow(beam_num_ptcl/math.pi/n_pe,3/8.)/pow(beam_tot_z,1/8.)
    return rb_max

def calc_power_beam_plasma(n_pe, rb_max):
    """
    Calculate the power transferred from the beam to the plasma bubble

    Valid in "strong bubble regime", where rb_max*k_pe >> 1.
    Args:
        n_pe:   number density of the electron plasma
        rb_max: maximum radius of the plasma bubble
    Returns:
        p_beam_plasma: power transferred from beam to plasma bubble
    """
    # from Eq. (10) of LBN2017
    p_beam_plasma = rsconst.c*(0.5*math.pi*n_pe*rsconst.e*rsconst.MKS_factor*rb_max**2)**2
    return p_beam_plasma


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
    # from Eq. (2) of LBN2017
    Ez = -2.*math.pi*n_pe*np.abs(rsconst.e*rsconst.MKS_factor)*rb*drb_dxi
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
    # from Eq. (3) of LBN2017
    # there is ambiguity in the sign, which needs to be resolved
    drb_dxi = math.sqrt((pow(rb_max/rb,4)-1.)/2.)
    return drb_dxi

def calc_Ez_on_axis_no_beam(n_pe, rb, rb_max):
    """
    Calculate the longitudinal electric field in a plasma bubble.

    Valid in "strong bubble regime", where rb_max*k_pe >> 1.
    Valid only in back of bubble, where there is no drive beam.
    The calculation is local; it depends on bubble radius and slope.
    Args:
        n_pe:    number density of the electron plasma
        rb:      the local bubble radius
        rb_max: maximum value of the bubble radius
    Returns:
        Ez: longitudinal electric field along the axis
    """
    drb_dxi = calc_drb_dxi_no_beam(rb, rb_max)

    # from Eq. (3) of LBN2017
    # there is ambiguity in the sign, which needs to be resolved
    Ez = -2.*math.pi*n_pe*np.abs(rsconst.e*rsconst.MKS_factor)*rb*drb_dxi
    return Ez

def calc_bubble_halfwidth(rb_max):
    """
    Calculate the halfwidth of the plasma bubble

    Valid in "strong bubble regime", where rb_max*k_pe >> 1.
    Args:
        rb_max: maximum value of the bubble radius
    Returns:
        rb: the local bubble radius
    """
    # from Eq. (4) of LBN2017
    xi_b = 0.847*rb_max
    return xi_b

def calc_local_bubble_radius(xi, rb_max):
    """
    Calculate the local bubble radius rb(xi)

    Valid in "strong bubble regime", where rb_max*k_pe >> 1.
    xi=ct-z, the distance from the front of the bubble (positive)
    Args:
        xi:     distance from front of the bubble
        rb_max: maximum value of the bubble radius
    Returns:
        rb: the local bubble radius
    """
    # halfwidth of the plasma bubble
    xi_b = calc_bubble_halfwidth(rb_max)

    # from Eq. (5) of LBN2017
    if xi>=2*xi_b or xi<=0.: rb = 0.
    else: rb = rb_max*math.pow((1.-((xi-xi_b)/xi_b)**2),(1./3.))
    return rb

def calc_E_decel_along_beam(n_pe, beam_tot_z, beam_num_ptcl):
    """
    Calculate (constant?) decelerating Ez along the beam (on axis)

    Valid in "strong bubble regime", where rb_max*k_pe >> 1.
    Args:
        n_pe:           number density of the electron plasma
        beam_tot_z:     total length of the drive beam
        beam_num_ptcl:  number of e- in the drive beam
    Returns:
        E_decel: (constant?) Ez along the beam (on axis)
    """
    # the following is large, when the calculation is valid
    strong_check_2 = beam_num_ptcl/n_pe/beam_tot_z**3

    # derived from Eq. (8) of LBN2017
    E_decel = math.pi * n_pe * beam_tot_z * \
              np.abs(rsconst.e * rsconst.MKS_factor) * \
              (math.sqrt(1. + 8. * strong_check_2 / math.pi) - 1.)
    return E_decel

