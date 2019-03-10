# -*- coding: utf-8 -*-
"""
Read openPMD field data from HDF5 file without using openPMD.

:copyright: Copyright (c) 2019 Radiasoft LLC. All Rights Reserved.
:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""

# Python imports
import h5py

def read_vector(path_to_file, field_name, field_coord, n_dump_str):
    """
    Read one component of a vector field from an HDF5 file.

    Assume 3 components (r,t,z) with openPMD conventions
    Assume 2D mesh of values (ie quasi-3D rz)
    Args:
        path_to_file: location of a specific HDF5 file
        field_name:   name of field in the HDF5 file
        field_coord:  field coordinate ('r','t', or 'z')
        n_dump_str:   dump number (as a string)
    Returns:
        field:     specified component of the field data
        sim_time:  time [s] at which data was dumped
        grid_size: grid spacings dr [m] and dz [m]
   """
    file = h5py.File(path_to_file,'r')
    data = file.get('data/')
    step = data.get(n_dump_str)
    all_fields = step.get('fields')
    my_field = all_fields.get(field_name)
    field_h5 = my_field.get(field_coord)
    field = field_h5[0,:,:]

    sim_time = step.attrs["time"] * step.attrs["timeUnitSI"]
    grid_size = my_field.attrs["gridSpacing"] * my_field.attrs["gridUnitSI"]

    return field, sim_time, grid_size

def read_scalar(path_to_file, field_name, n_dump_str):
    """
    Read a scalar field by name from an HDF5 file.

    Assume openPMD conventions
    Assume 2D mesh of values (ie quasi-3D rz)
    Args:
        path_to_file: location of a specific HDF5 file
        field_name:   name of field in the HDF5 file
        n_dump_str:   dump number (as a string)
    Returns:
        field:     the requested field data
        sim_time:  time [s] at which data was dumped
        grid_size: grid spacings dr [m] and dz [m]
   """
    file = h5py.File(path_to_file,'r')
    data = file.get('data/')
    step = data.get(n_dump_str)
    all_fields = step.get('fields')
    field_h5 = all_fields.get(field_name)
    field = field_h5[0,:,:]

    sim_time = step.attrs["time"] * step.attrs["timeUnitSI"]
    grid_size = field_h5.attrs["gridSpacing"] * field_h5.attrs["gridUnitSI"]

    return field, sim_time, grid_size
