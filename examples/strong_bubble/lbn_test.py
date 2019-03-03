from __future__ import absolute_import, division, print_function, unicode_literals

import math
import numpy as np

from rsbeams.rsstats.stats6d import specify_significant_figures as rs_sigfig
from rsbeams.rsphysics import rsconst
from rsfbpic.rswake import lbn_wake

# Specify default values of some physical quantities...
#   These are motivated by upcoming PWFA experiments at FACET II

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
# beam_rms_r = 3.65e-6   # [m]   FACET II params...?
beam_rms_r = 0.4 / k_pe  # [m]   satisfies resonance condition

# RMS length
# beam_rms_z = 12.77e-6  # [m]   FACET II params...?
beam_rms_z = 1.2 / k_pe  # [m]   satisfies resonance condition

# total charge in the electron drive beam
# beam_tot_q = 1.e10*np.abs(rsconst.e)    # [C]   FACET II params...?
beam_tot_q = 3.e-9                          # [C]   3 nC
# peak number density of the drive beam
beam_rms_v = math.pi*beam_rms_z*beam_rms_r**2
beam_rms_n = beam_tot_q/np.abs(rsconst.e)/beam_rms_v

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

# -----------------
# Excercise some of the methods
#------------------

rb = 0.4 * lambda_pe
drb_dxi = 2.73861278753
Ez = lbn_wake.calc_Ez_on_axis(n_pe, rb, drb_dxi)
print()
print("Ez = ", rs_sigfig(Ez*1.e-9,3), " [GV/m]")

rb_max = 2.*rb
drb_dxi = lbn_wake.calc_drb_dxi_no_beam(rb, rb_max)
print()
print("drb_dxi = ", rs_sigfig(drb_dxi,3), " [rad]")


print()
