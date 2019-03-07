# -*- coding: utf-8 -*-
"""
Read openPMD field data from HDF5 file without using openPMD.

:copyright: Copyright (c) 2019 Radiasoft LLC. All Rights Reserved.
:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""

# Python imports
import h5py

def read_vector(path_to_file, field_name, f_coord):
    """
    Read in a field by name from an HDF5 file.

    Args:
        path_to_file: location of a specific HDF5 file
        field_name:   name of field in the HDF5 file
    Returns:
        data:  TBD...
    """
    file = h5py.File(path_to_file,'r')
    data = file.get('data/')
    step = data.get('280')
    fields = step.get('fields')
    Es = fields.get('E')
    Ez = Es.get('z')
    theEz = Ez[0,:,:]

    return theEz
