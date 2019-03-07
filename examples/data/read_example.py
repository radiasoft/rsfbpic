# Python imports
from __future__ import absolute_import, division, print_function, unicode_literals

# RadiaSoft imports
from rsfbpic.rsopenpmd import read_vector_field

# -------------------
# read some data
# -------------------

# read in the longitudinal electric field
path_to_data = '../../rsfbpic/package_data'
field_name = 'E'
f_coord = 'z'
n_dump = 280
mode_number = 0
theta_value = 0.
show_plot = False

ez, info_ez = read_vector_field.read_vector(path_to_data, field_name, f_coord, n_dump, mode_number, theta_value, show_plot)

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

print()
print("Ez[r=0,z] = ", ez[64,:])
