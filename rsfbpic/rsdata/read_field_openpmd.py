# -*- coding: utf-8 -*-
"""
Read field data from HDF5 file via OpenPMD standard.

:copyright: Copyright (c) 2019 Radiasoft LLC. All Rights Reserved.
:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""

# OpenPMD imports
from opmd_viewer import OpenPMDTimeSeries

def read_vector(path_to_data, field_name, f_coord, n_dump, \
                mode_number, theta_value, show_plot):
    """
    Read in a field by name from an HDF5 file.

    Args:
        path_to_data: location of directory for HDF5 files
        field_name:   name of field in the HDF5 files
        n_dump:       dump number to be read
    Returns:
        field:  gridded r-z field data
        info_f: metadata regarding the field
    """

    # open the file; instantiate "time series" object
    time_series = OpenPMDTimeSeries(path_to_data)

    # read the specified field
    field, info_f = time_series.get_field(iteration=n_dump, \
                                          field=field_name, \
                                          coord=f_coord,    \
                                          m=mode_number,    \
                                          theta=theta_value,\
                                          plot=show_plot)
    return field, info_f
