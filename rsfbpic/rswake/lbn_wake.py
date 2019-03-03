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


def calc_bubble_half_width(rb_max):
    """
    Calculate the half width of the plasma bubble
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
    # half width of the plasma bubble
    xi_b = calc_bubble_half_width(rb_max)
    

    # from Eq. (5) of LBN2017
    if xi>=2*xi_b or xi<=0.: rb = 0.
    else: rb = rb_max*math.pow((1.-((xi-xi_b)/xi_b)**2),(1./3.))
    return rb
