"""
This is an input script for modeling a PWFA with an axi-symmetric drive beam,
with or without ion motion included, and on many cores.
"""

# Imports
# standard python libraries
import numpy as np
from scipy import constants
from scipy.special import erfc, k0, k1

import shutil, os

import matplotlib.pyplot as plt
import matplotlib as mpl

import h5py as hdf5

# Imports for the simulations, and setting up the plots
from fbpic.main import Simulation
from fbpic.openpmd_diag import FieldDiagnostic, ParticleDiagnostic, \
     set_periodic_checkpoint, restart_from_checkpoint # ParticleChargeDensityDiagnostic,\
from fbpic.lpa_utils.bunch import add_elec_bunch_gaussian

# set the colormap and centre the colorbar

import matplotlib.colors as colors

# function for centering the colorbars so the middle color always corresponds to 0.
class MidpointNormalize(colors.Normalize):
    """
    Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

####################
##
## Physical parameters
##
####################


## Beam parameters

# Drive bunch is gaussian

# turn the drive beam on and off, as required.
use_drive_bunch = True

drive_sigma_r = 3.65e-6  # meters
drive_sigma_z = 12.77e-6  # meters
drive_Q = 1.e10*(-1.*constants.elementary_charge)   # Coulombs
drive_N_macro = 10000
# We use a very large gamma to avoid having to deal with an evolving
# drive bunch, which makes computing wake functions more complicated
drive_gamma = 10.e9

# Witness bunch, also gaussian

# turn the witness beam on and off, as required.
use_witness_bunch = True

witness_sigma_r = 3.65e-6 #meters
witness_sigma_z = 6.38e-6  # meters
witness_Q = 4.3e9*(-1.*constants.elementary_charge)   # Coulombs
witness_N_macro = 7500
# We use a very large gamma to avoid having to deal with an evolving
# witness bunch, which makes computing wake functions more complicated
witness_gamma = 10.e9
trailing_distance = 150.e-6 # meters

# flag that includes ion motion
# Set to true, this will create a charge-neutral plasma with equal densities
# of plasma electrons and singly ionized ions.
include_ion_motion = True

# for now, assume Helium
ion_mass = 2.*constants.proton_mass + 2.*constants.neutron_mass
ion_charge = constants.elementary_charge


## Plasma channel parameters

n_plasma = 4.e16        # cm^-3

# convert to per cubic meter
n_plasma *= 100**3

# derived plasma quantities
omega_p = np.sqrt(n_plasma*constants.elementary_charge**2/(constants.m_e*constants.epsilon_0))
k_p = omega_p/constants.c

lambda_p = 2.*np.pi/k_p

####################
##
## Simulation parameters
##
####################

## Number of macroparticles per cell

n_macro_r = 4
n_macro_z = 4
n_macro_t = 1 # this should be 1 if we are not using higher order azimuthal modes

# where to dump the data
dump_dir = './diags'

# We want to run the simulation just long enough for the fields to form behind the drive bunch,

# Domain size, include the whole initial bucket and some trailing distance
domain_length = 3.*lambda_p  # meters
domain_radius = lambda_p  # meters

# start the ramp after the drive bunch has existed a while
ramp_start = domain_length
ramp_length = 5.*drive_sigma_z

# Grid size, resolve the drive bunch
Delta_z = min([0.05*drive_sigma_z, 0.05*lambda_p])  # meters
Delta_r = min([0.05*drive_sigma_r, 0.05*lambda_p])  # meters

# Derived quantities
Nz = int(np.rint(domain_length/Delta_z))
Nr = int(np.rint(domain_radius/Delta_r))

# One cell per step, for the moving window
dt = (domain_length)/constants.c/Nz  # sec

sim_length = (ramp_start + ramp_length + 3*domain_length)/constants.c

Nsteps = int(sim_length/dt)-int(sim_length/dt)%100 + 1

## Define the diagnostics
write_fields = True
write_particles = True
# In this context, we only look at the fields at the final step
dump_period = Nsteps-1

# Moving window
window_v = constants.c

# n_order specifies the deposition stencil. -1 is global and works in serial, and gets
# exactly correct dispersion, but cannot run in parallel. A rule of thumb from the
# fbpic user community is to use 32-cell longitudinal deposition as a trade-off between
# fidelity and
n_order = -1

# create the density function for the plasma, which is uniform
def dens_func( z, r ) :
    """Returns relative density at position z and r"""
    # Allocate relative density
    n = np.ones_like(z)
    # Make linear ramp
    n = np.where( z < ramp_start + ramp_length, (z-ramp_start)/ramp_length, n )
    # Supress density before the ramp
    n = np.where( z < domain_length + ramp_start, 0., n )
    return(n)

# Use only the primary azimuthal mode.
Nm = 1

####################
##
## Create and run the simulation
##
## This section shouldn't be changed from simulation to simulation
##
####################

# remove old data
# This means you should probably change dump_dir if you are running a battery of simulations
if os.path.exists(dump_dir):
    shutil.rmtree(dump_dir)

# Create the simulation
sim = Simulation(Nz, domain_length, Nr, domain_radius, Nm, dt,
                 boundaries='open', n_order=n_order)
# micromanage the particle species by removing the default created species first
sim.ptcl = []

# add the gaussian drive bunch
if use_drive_bunch:
    add_elec_bunch_gaussian( sim,
                            sig_r = drive_sigma_r,
                            sig_z = drive_sigma_z,
                            n_emit=0.,
                            gamma0=drive_gamma,
                            sig_gamma=1.,
                            Q=drive_Q,
                            N=drive_N_macro,
                            tf=0.0,
                            zf=.75*domain_length, boost=None)
if use_witness_bunch:
    add_elec_bunch_gaussian( sim,
                            sig_r = witness_sigma_r,
                            sig_z = witness_sigma_z,
                            n_emit=0.,
                            gamma0=witness_gamma,
                            sig_gamma=1.,
                            Q=witness_Q,
                            N=witness_N_macro,
                            tf=0.0,
                            zf=.75*domain_length - trailing_distance, boost=None)

# add the plasma electrons
plasma_electrons = sim.add_new_species(q = -1.*constants.elementary_charge,
                                 m = constants.electron_mass,
                                 dens_func = dens_func,
                                 n = n_plasma,
                                 p_nz = n_macro_z, p_nr = n_macro_r, p_nt = n_macro_t)
if include_ion_motion:
    plasma_ions = sim.add_new_species(q = ion_charge,
                                 m = ion_mass,
                                 dens_func = dens_func,
                                 n = n_plasma,
                                 p_nz = n_macro_z, p_nr = n_macro_r, p_nt = n_macro_t)

# Set the moving window
sim.set_moving_window(v = window_v)

# Add diagnostics
if write_fields:
    sim.diags.append( FieldDiagnostic(dump_period, sim.fld, sim.comm, write_dir=dump_dir ) )
if write_particles:
    sim.diags.append( ParticleDiagnostic( dump_period,
                    {'electrons': sim.ptcl[0]}, sim.comm, write_dir=dump_dir ) )

# run the simulation
sim.step(Nsteps)
