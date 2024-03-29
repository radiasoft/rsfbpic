# Python imports
from __future__ import absolute_import, division, print_function, unicode_literals

# RadiaSoft imports
from rsfbpic.rsdata import read_field_hdf

# -------------------
# read data via h5py
# -------------------

print()
print("---- read E fields -----")
print()

# read in the longitudinal electric field component
path_to_file = '../../rsfbpic/package_data/data00000280.h5'
field_name = 'E'
field_coord = 'z'
n_dump = 280
n_dump_str = str(n_dump)

time = read_field_hdf.read_time(path_to_file, n_dump_str)

dr, dz = read_field_hdf.read_dr_dz(path_to_file, field_name, n_dump_str)

ez = read_field_hdf.read_vector(path_to_file, field_name, field_coord, n_dump_str)

print()
print("simulation time = ", time, " [s]")
print("radial grid sized = ", dr, " [m]")
print("axial grid sized = ", dz, " [m]")

print()
print()

# normalize units to GV/m and microns
ez *= 1.e-9
# xi_array *= 1.e6

# generate the plot
# ax = plt.subplot(111)
# ax.plot(xi_array, ez)
# ax.set_xlabel('xi = ct - z [microns]')
# ax.set_ylabel('(axial) Ez [GV/m]')
# ax.set_title('PWFA axial Ez in "strong" regime')
