from __future__ import absolute_import, division, print_function, unicode_literals
import pytest

import math
import numpy as np

from rsbeams.rsstats.stats6d import specify_significant_figures as rs_sigfig
from rsbeams.rsphysics import rsconst
from rsfbpic.rswake import lbn_wake

# Specify values of relevant physical quantities...

# number density of pre-ionized electron plasma
n_pe_cgs = 4.e16             # [cm^-3]
n_pe     = n_pe_cgs * 1.e6   #  [m^-3]

# electron plasma frequency, wavenumber, wavelength
om_pe     = np.sqrt(n_pe*rsconst.e**2
            / (rsconst.m_e*rsconst.epsilon_0) )  # [rad/s]
k_pe      = om_pe/rsconst.c                      # [rad/m]
lambda_pe = 2.*np.pi/k_pe                        # [m]
lambda_pe_microns = lambda_pe*1.e6               # [microns]

# -------------------
# The PWFA drive beam, with an assumed Gaussian distribution
# -------------------

# RMS radius
beam_rms_r = 0.4 / k_pe    # [m] satisfies resonance condition
beam_max_r = 3*beam_rms_r  # [m] ad hoc definition of max radius

# RMS length
beam_rms_z = 0.5 / k_pe    # [m] satisfies resonance condition
beam_tot_z = 4*beam_rms_z  # [m] ad hoc definition of total length

# total charge in the electron drive beam
beam_tot_q = 3.e-9                           # [C] 3 nC
beam_num_ptcl = beam_tot_q/np.abs(rsconst.e) # [] # of beam electrons

# peak number density of the drive beam
beam_rms_v = math.pi*beam_rms_z*beam_rms_r**2
beam_rms_n = beam_tot_q/np.abs(rsconst.e)/beam_rms_v

# large density ratio is required for blowout regime
dens_ratio = beam_rms_n/n_pe    # should be >> unity

# relativistic gamma factor
beam_gamma = 1.957e+4    # FACET II params...?

# maximum radius of the plasma bubble
def test_lbn_01():
    rb_max = lbn_wake.calc_rb_max(n_pe, beam_tot_z, beam_num_ptcl)
    assert rs_sigfig(rb_max*1.e6,3) == rs_sigfig(97.2,3)

# the 1st 'strong bubble' validity condition
def test_lbn_02():
    rb_max = lbn_wake.calc_rb_max(n_pe, beam_tot_z, beam_num_ptcl)
    strong_check_1 = rb_max*k_pe
    assert rs_sigfig(strong_check_1,3) == rs_sigfig(3.66,3)

# the 2nd 'strong bubble' validity condition
def test_lbn_03():
    strong_check_2 = beam_num_ptcl/n_pe/beam_tot_z**3
    assert rs_sigfig(strong_check_2,3) == rs_sigfig(3.12,3)

# power transferred from the beam to the plasma
def test_lbn_04():
    rb_max = lbn_wake.calc_rb_max(n_pe, beam_tot_z, beam_num_ptcl)
    p_beam_plasma = lbn_wake.calc_power_beam_plasma(n_pe, rb_max)
    assert rs_sigfig(p_beam_plasma,3) == rs_sigfig(2.19e+20,3)

# local axial longitudinal electric field
def test_lbn_05():
    rb = 0.4 * lambda_pe
    drb_dxi = 2.73861278753
    Ez = lbn_wake.calc_Ez_on_axis(n_pe, rb, drb_dxi)
    assert rs_sigfig(Ez*1.e-9,3) == rs_sigfig(-66.2,3)

# local derivative of the bubble radius
def test_lbn_06():
    rb = 0.4 * lambda_pe
    rb_max = lbn_wake.calc_rb_max(n_pe, beam_tot_z, beam_num_ptcl)
    drb_dxi_nb = lbn_wake.calc_drb_dxi_no_beam(rb, rb_max)
    assert rs_sigfig(drb_dxi_nb,3) == rs_sigfig(1.32,3)

# local axial longitudinal electric field
def test_lbn_07():
    rb = 0.4 * lambda_pe
    rb_max = lbn_wake.calc_rb_max(n_pe, beam_tot_z, beam_num_ptcl)
    Ez_nb = lbn_wake.calc_Ez_on_axis_no_beam(n_pe, rb, rb_max)
    assert rs_sigfig(Ez_nb*1.e-9,3) == rs_sigfig(-31.9,3)

# halfwidth of the plasma bubble
def test_lbn_08():
    rb_max = lbn_wake.calc_rb_max(n_pe, beam_tot_z, beam_num_ptcl)
    xi_b = lbn_wake.calc_bubble_halfwidth(rb_max)
    assert rs_sigfig(xi_b*1.e6,3) == rs_sigfig(82.3,3)

# local bubble radius at arbitrary location
def test_lbn_09():
    rb_max = lbn_wake.calc_rb_max(n_pe, beam_tot_z, beam_num_ptcl)
    xi_b = lbn_wake.calc_bubble_halfwidth(rb_max)
    xi = 0.3*xi_b
    rb = lbn_wake.calc_local_bubble_radius(xi, rb_max)
    assert rs_sigfig(rb*1.e6,3) == rs_sigfig(77.7,3)

# bubble radius at location of maximum
def test_lbn_10():
    rb_max = lbn_wake.calc_rb_max(n_pe, beam_tot_z, beam_num_ptcl)
    xi_b = lbn_wake.calc_bubble_halfwidth(rb_max)
    rb = lbn_wake.calc_local_bubble_radius(xi_b, rb_max)
    assert rs_sigfig(rb*1.e6,3) == rs_sigfig(97.2,3)

# decelerating E field along the axis
def test_lbn_11():
    E_decel = lbn_wake.calc_E_decel_along_beam(n_pe, beam_tot_z, beam_num_ptcl)
    assert rs_sigfig(E_decel*1.e-9,3) == rs_sigfig(19.1,3)
