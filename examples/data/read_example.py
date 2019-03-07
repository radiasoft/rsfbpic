# Python imports
from __future__ import absolute_import, division, print_function, unicode_literals

# RadiaSoft imports
from rsfbpic.rsdata import read_field_hdf
from rsfbpic.rsdata import read_field_openpmd

# -------------------
# read data via openPMD
# -------------------

# read in the longitudinal electric field
path_to_data = '../../rsfbpic/package_data'
field_name = 'E'
f_coord = 'z'
n_dump = 280
mode_number = 0
theta_value = 0.
show_plot = False

ez, info_ez = read_field_openpmd.read_vector(path_to_data, field_name, f_coord, n_dump, mode_number, theta_value, show_plot)

print()
print("---- using openPMD -----")
print()
print("axes = ", info_ez.axes)
print()
print("rmin, rmax = ", info_ez.rmin, "; ", info_ez.rmax)
print("zmin, zmax = ", info_ez.zmin, "; ", info_ez.zmax)
print()
print("dz = ", info_ez.dz)
print("dr = ", info_ez.dr)
print()
# print("r = ", info_ez.r)
# print("z = ", info_ez.z)
print("N_cells_r = ", len(info_ez.r))
print("N_cells_z = ", len(info_ez.z))

# print()
# print("Ez[r=0,z] = ", ez[64,:])


# -------------------
# read data via h5py
# -------------------

print()
print("---- using h5py -----")
print()

# read in the longitudinal electric field
path_to_file = '../../rsfbpic/package_data/data00000280.h5'
field_name = 'E'
f_coord = 'z'
n_dump = 280

ez2 = read_field_hdf.read_vector(path_to_file, field_name, f_coord)

print("length of ez2 = ", len(ez2))
print("ez2 = ", ez2)

print()
print()

# normalize units to GV/m and microns
ez_array *= 1.e-9
xi_array *= 1.e6

# generate the plot
ax = plt.subplot(111)
ax.plot(xi_array, ez_array)
ax.set_xlabel('xi = ct - z [microns]')
ax.set_ylabel('(axial) Ez [GV/m]')
ax.set_title('PWFA axial Ez in "strong" regime')
