from rsbeams.rsdata.SDDS import writeSDDS
import numpy as np

# Create some example data. Here we mimic creating the 6D particle distribution of an electron bunch
# in a particle accelerator

coordinates = np.random.rand(100, 6)
total_particles = 100
central_momentum = 0.984

# Intialize the SDDS writer
new_distr = writeSDDS()

# Add new columns - Name, data, and data type are required. See the `create_column` docstring for optional parameters
new_distr.create_column('x', coordinates[:, 0], 'double', colUnits='m')
new_distr.create_column('xp', coordinates[:, 1], 'double', colUnits='')
new_distr.create_column('y', coordinates[:, 2], 'double', colUnits='m')
new_distr.create_column('yp', coordinates[:, 3], 'double', colUnits='')
new_distr.create_column('t', coordinates[:, 4], 'double', colUnits='s')
new_distr.create_column('p', coordinates[:, 5], 'double', colUnits='')

# Add new parameters - Name, data, and data type are required.
# See the `create_parameter` docstring for optional parameters
new_distr.create_parameter('particles', total_particles, 'long')
new_distr.create_parameter('pCentral', central_momentum, 'double', parUnits='MeV/c')

# Save the file. SDDS binary and ascii formats are supported.
new_distr.save_sdds('example_sdds_bunch.sdds', dataMode='binary')