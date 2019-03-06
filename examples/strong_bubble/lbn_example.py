from __future__ import absolute_import, division, print_function, unicode_literals

import math
import numpy as np
import matplotlib.pyplot as plt

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
# beam_rms_r = 3.65e-6     # [m] FACET II params...?
beam_rms_r = 0.4 / k_pe    # [m] satisfies resonance condition
beam_max_r = 3*beam_rms_r  # [m] ad hoc definition of max radius

# RMS length
# beam_rms_z = 12.77e-6    # [m] FACET II params...?
beam_rms_z = 0.5 / k_pe    # [m] satisfies resonance condition
beam_tot_z = 4*beam_rms_z  # [m] ad hoc definition of total length

# total charge in the electron drive beam
# beam_tot_q = 1.e10*np.abs(rsconst.e)       # [C] FACET II params...?
beam_tot_q = 3.e-9                           # [C] 3 nC
beam_num_ptcl = beam_tot_q/np.abs(rsconst.e) # [] # of beam electrons

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

print()
print("*******")
print("Calculate the maximum radius of the plasma bubble:")
rb_max = lbn_wake.calc_rb_max(n_pe, beam_tot_z, beam_num_ptcl)
print("    rb_max = ", rs_sigfig(rb_max*1.e6,3), " [microns]")

print()
print("*******")
print("Calculate the 1st 'strong bubble' validity condition...")
print("(bubble radius) / (plasma skin depth); it must be large:")
strong_check_1 = rb_max*k_pe
print("    1st validity ratio = ", rs_sigfig(strong_check_1,3))

print()
print("*******")
print("Calculate the 2nd 'strong bubble' validity condition...")
print("(scaled beam dens) / (plasma dens); it must be large:")
strong_check_2 = beam_num_ptcl/n_pe/beam_tot_z**3
print("    2nd validity ratio = ", rs_sigfig(strong_check_2,3))

print()
print("*******")
print("Calculate power transferred from the beam to the plasma:")
p_beam_plasma = lbn_wake.calc_power_beam_plasma(n_pe, rb_max)
print("    p_beam_plasma = ", rs_sigfig(p_beam_plasma,3), " [W]")

print()
print("*******")
print("Calculate the local axial longitudinal electric field:")
rb = 0.4 * lambda_pe
drb_dxi = 2.73861278753
Ez = lbn_wake.calc_Ez_on_axis(n_pe, rb, drb_dxi)
print("    Ez = ", rs_sigfig(Ez*1.e-9,3), " [GV/m]")

print()
print("*******")
print("Calculate the local derivative of the bubble radius...")
print("    We assume this is a location behind the drive beam:")
# rb_max is calculated above
drb_dxi_nb = lbn_wake.calc_drb_dxi_no_beam(rb, rb_max)
print("    drb_dxi (no beam) = ", rs_sigfig(drb_dxi_nb,3), " [rad]")

print()
print("*******")
print("Calculate local axial longitudinal electric field...")
print("    We assume this is a location behind the drive beam:")
Ez_nb = lbn_wake.calc_Ez_on_axis_no_beam(n_pe, rb, rb_max)
print("    Ez (no beam) = ", rs_sigfig(Ez_nb*1.e-9,3), " [GV/m] (+ or -)")

print()
print("*******")
print("Calculate the halfwidth of the plasma bubble:")
xi_b = lbn_wake.calc_bubble_halfwidth(rb_max)
print("    xi_b = ", rs_sigfig(xi_b*1.e6,3), " [microns]")

print()
print("*******")
print("Calculate local bubble radius at arbitrary location...")
print("    Selected location is 'xi=ct-z' (positive):")
xi = 0.3*xi_b
print("    xi = ", rs_sigfig(xi*1.e6,3), " [microns]")
rb = lbn_wake.calc_local_bubble_radius(xi, rb_max)
print("    rb = ", rs_sigfig(rb*1.e6,3), " [microns]")

print()
print("*******")
print("Calculate bubble radius at location of maximum...")
print("    Selected location is 'xi=xi_b' (positive):")
rb = lbn_wake.calc_local_bubble_radius(xi_b, rb_max)
print("    rb = ", rs_sigfig(rb*1.e6,3), " [microns]")
print("    rb/rb_max = ", rs_sigfig(rb/rb_max,3), " (should be unity!)")

print()
print("*******")
print("Calculate the decelerating E field along the axis...")
print("    Assumed constant ??, it is along the drive beam:")
E_decel = lbn_wake.calc_E_decel_along_beam(n_pe, beam_tot_z, beam_num_ptcl)
print("    E_decel = ", rs_sigfig(E_decel*1.e-9,3), " [GV/m]")

print()
print("*******")

# now we'll create some plots

from rsfbpic.rsplot import lbn_plot
import matplotlib.pyplot as plt

# plot the axial electric field
ez_plot = lbn_plot.plot_Ez_on_axis(n_pe, beam_tot_z, beam_num_ptcl)
plt.show()

# plot the bubble radius
rb_plot = lbn_plot.plot_bubble_radius(n_pe, beam_tot_z, beam_num_ptcl)
plt.show()
