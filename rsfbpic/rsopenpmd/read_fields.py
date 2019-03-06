# -*- coding: utf-8 -*-
"""
Plotting curves from 2017 PRAB article by Lebedev, Burov and Nagaitsev

:copyright: Copyright (c) 2019 Radiasoft LLC. All Rights Reserved.
:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""

# SciPy imports
import numpy as np
import matplotlib.pyplot as plt

# RadiaSoft imports
from rsfbpic.rswake import lbn_wake

def plot_Ez_on_axis(n_pe, beam_tot_z, beam_num_ptcl):
    """
    Plot the longitudinal electric field in a plasma bubble.

    Valid in "strong bubble regime", where rb_max*k_pe >> 1.
    Args:
        n_pe:    number density of the electron plasma
        beam_tot_z:     total length of the drive beam
        beam_num_ptcl:  number of e- in the drive beam
    Returns:
        ax: matplotlib 'Axis' object for Ez inside bubble
    """

    # calculate the maximum bubble radius
    rb_max = lbn_wake.calc_rb_max(n_pe, beam_tot_z, beam_num_ptcl)

    # calculate the bubble half width
    xi_b = lbn_wake.calc_bubble_halfwidth(rb_max)

    # calculate the approximate (constant?) decelerating field
    #     along the drive beam
    E_decel = lbn_wake.calc_E_decel_along_beam(n_pe, beam_tot_z, beam_num_ptcl)

    # Specify the plot range, xi_min <= xi <= xi_max
    # xi=ct-z is the distance from the front of the bubble (positive)
    xi_min = 0.
    xi_max = 1.99*xi_b
    num_points = 100
    xi_array = np.linspace(xi_min, xi_max, num=num_points)
    ez_array = np.zeros(num_points)
    for iloop in range(0, num_points):
        xi = xi_array[iloop]
        if xi < 0. or xi > 2.*xi_b: ez_array[iloop] = 0.
        elif xi < beam_tot_z: ez_array[iloop] = E_decel
        else:
            rb = lbn_wake.calc_local_bubble_radius(xi, rb_max)
            ez_array[iloop] = lbn_wake.calc_Ez_on_axis_no_beam(n_pe, rb, rb_max)

    # normalize units to GV/m and microns
    ez_array *= 1.e-9
    xi_array *= 1.e6

    # generate the plot
    ax = plt.subplot(111)
    ax.plot(xi_array, ez_array)
    ax.set_xlabel('xi = ct - z [microns]')
    ax.set_ylabel('(axial) Ez [GV/m]')
    ax.set_title('PWFA axial Ez in "strong" regime')
    return ax

def plot_bubble_radius(n_pe, beam_tot_z, beam_num_ptcl):
    """
    Plot the plasma bubble radius.

    Valid in "strong bubble regime", where rb_max*k_pe >> 1.
    Args:
        n_pe:    number density of the electron plasma
        beam_tot_z:     total length of the drive beam
        beam_num_ptcl:  number of e- in the drive beam
    Returns:
        ax: matplotlib 'Axis' object for bubble radius plot
    """

    # calculate the maximum bubble radius
    rb_max = lbn_wake.calc_rb_max(n_pe, beam_tot_z, beam_num_ptcl)

    # calculate the bubble half width
    xi_b = lbn_wake.calc_bubble_halfwidth(rb_max)

    # Specify the plot range, xi_min <= xi <= xi_max
    # xi=ct-z is the distance from the front of the bubble (positive)
    xi_min = 0.
    xi_max = 2.*xi_b
    num_points = 200
    xi_array = np.linspace(xi_min, xi_max, num=num_points)
    rb_array = np.zeros(num_points)
    for iloop in range(0, num_points):
        rb_array[iloop] = lbn_wake.calc_local_bubble_radius(xi_array[iloop], rb_max)

    # normalize units to microns
    rb_array *= 1.e6
    xi_array *= 1.e6

    # generate the plot
    ax = plt.subplot(111)
    ax.plot(xi_array, rb_array)
    ax.set_xlabel('xi = ct - z [microns]')
    ax.set_ylabel('rb [microns]')
    ax.set_title('PWFA bubble radius in "strong" regime')
    return ax
