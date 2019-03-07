# Python imports
from __future__ import absolute_import, division, print_function, unicode_literals

# RadiaSoft imports
from rsfbpic.rsopenpmd import read_vector_field

# -------------------
# read some data
# -------------------

# read in the electric field
path_to_data = '../../rsfbpic/package_data'
field_name = 'E'
f_coord = 'z'
n_dump = 280
e_fields = read_vector_field.read_vector(path_to_data, field_name, f_coord, n_dump)

