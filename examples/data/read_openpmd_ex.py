# Python imports
from __future__ import absolute_import, division, print_function, unicode_literals

# OpenPMD imports
from opmd_viewer import OpenPMDTimeSeries

# read in the longitdunal electric field component
path_to_data = '../../rsfbpic/package_data'
field_name = 'E'
coord = 'z'
iteration = 280
mode = 0
theta = 0.
plot = False

# open the file; instantiate "time series" object
time_series = OpenPMDTimeSeries(path_to_data)

# read the specified field
ez, info_ez = time_series.get_field(iteration=iteration, \
                                    field=field_name, \
                                    coord=coord, \
                                    m=mode, \
                                    theta=theta, \
                                    plot=plot)
                
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
