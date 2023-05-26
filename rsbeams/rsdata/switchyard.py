
# Import the relevant data formats

import pandas as pd
import numpy as np
import h5py as h5
from scipy import constants
from rsbeams.rsptcls.species import Species
from subprocess import Popen, PIPE
from rsbeams.rsdata.SDDS import writeSDDS, readSDDS

_DEFAULT_SPECIES_NAME = 'species_0'


def read_elegant(file_name, species_name ='Species'):
    """Read in a file from elegant output.
    :file_name: name of file to read from
    """

    # elegant coordinates are as follows:
    # x  -- horizontal offset, m
    # x' -- horizontal angle, rad
    # y  -- vertical offset, m
    # y' -- vertical angle, rad
    # t  -- time of flight, sec
    # p  -- longitudinal momentum, mc

    read_particle_data = readSDDS(file_name)
    read_particle_data.read()
    sdds_data = read_particle_data.columns.squeeze()
    try:
        charge_data = read_particle_data.parameters['Charge'].squeeze()
    except ValueError:
        print(f"No Charge in file {file_name}")
        charge_data = 0.0

    particle_data = np.zeros([sdds_data.size, 6])
    for i, col in enumerate(['x', 'xp', 'y', 'yp', 't', 'p']):
        particle_data[:, i] = sdds_data[col][:]

    # Coordinate conversions
    particle_data[:, 1] *= particle_data[:, 5]
    particle_data[:, 3] *= particle_data[:, 5]
    particle_data[:, 4] *= constants.c

    species = Species(particle_data, charge=-1, mass=0.511e6, total_charge=charge_data)


    return species


def read_opal(file_name, step_number=None, species_name='Species'):
    """Read in a file from OPAL output.
    :file_name: name of file to read from
    """

    # opal coordinates are as follows:
    # x  -- horizontal offset, m
    # xp -- horizontal momentum, beta*gamma
    # y  -- vertical offset, m
    # yp -- vertical momentum, beta*gamma
    # z  -- Position relative to ? (some sort of reference), m
    # p  -- total momentum, beta*gamma

    with h5.File(file_name, 'r') as pcdata:
        if not step_number:
            # Use last step in the file
            all_step_numbers = []
            for key in pcdata.keys():
                sn_string = _get_step_number(key)
                try:
                    sn = int(sn_string)
                    all_step_numbers.append(sn)
                except ValueError:
                    pass
            step_number = np.max(all_step_numbers)
        loc = 'Step#{}'.format(step_number)
        mp_count = pcdata[loc+'/z'].shape[0]
        particle_data = np.empty((mp_count, 6))
        for i, coord in enumerate(['x', 'px', 'y', 'py', 'z', 'pz']):
            particle_data[:, i] = pcdata[loc+'/'+coord]
        total_charge = pcdata[loc].attrs['CHARGE'][0]


    # TODO: This shouldn't be specific to electrons
    species =  Species(particle_data, charge=-1, mass=0.511e6, total_charge=total_charge)

    return species


def read_genesis(file_name):
    """Read in a file from genesis output.
    :file_name: name of file to read from
    """

    # genesis coordinates are as follows:
    # x  -- horizontal coordinate, m
    # px -- horizontal momentum, beta_x gamma
    # y  -- vertical coordinate, m
    # py -- vertical momentum, beta_y gamma
    # th -- theta, the ponderomotive phase at fixed z
    # t  -- time of arrival at fixed z, sec
    # gamma -- particle gamma

    return 0


class Switchyard:
    """Class for writing a particle species data from a code output to a universal format, or vice verse. 
    Can be used to load output from one code into another, or for universal data visualization."""
    _supported_codes = ['genesis', 'elegant', 'opal']
    def __init__(self):
        
        # conventions
        # x  -- horizontal position, m
        # ux -- horizontal velocity, beta_x gamma
        # y  -- vertical position, m
        # uy -- vertical velocity, beta_y gamma
        # ct -- time of flight, m
        # pt -- total momentum, beta gamma
        self.species = {}
        self._readers = {'elegant': read_elegant, 'opal': read_opal}
        self._writers = {'elegant': self.write_elegant, 'genesis': self.write_genesis, 'opal': self.write_opal}
    
    def is_supported(self, code_name):
        
        if code_name in self.supported_codes():
            return True
        else:
            return False

    def _species_insert(self, species_name, species_object):
        result = self.species.setdefault(species_name, species_object)
        if result != species_object:
            print(f"Species named f{species_name} already exists.")

    @classmethod
    def supported_codes(cls):
        return cls._supported_codes

    def _get_reader(self, file_format):
        reader = self._readers.get(file_format)
        if not reader:
            supported = ', '.join(self.supported_codes())
            raise LookupError(f'Format {file_format} is not recognized. Available formats are: {supported}')

        return reader


    def read(self, file_name, file_format, species_name=None, **kwargs):
        reader = self._get_reader(file_format)
        if not species_name:
            species_name = _DEFAULT_SPECIES_NAME
        species = reader(file_name, species_name=species_name, **kwargs)
        self._species_insert(species_name, species)

    def write_elegant(self, file_name, species_name):
        """Write a file to elegant-readable format.
        :file_name: name of file to write to
        """
        
        # TODO Split out conversion function so we have converters and writers
        x = self.species[species_name].x
        y = self.species[species_name].y
        xp = self.species[species_name].ux / self.species[species_name].pt
        yp = self.species[species_name].uy / self.species[species_name].pt
        t = self.species[species_name].ct / constants.c
        p = self.species[species_name].pt
        
        file_out = writeSDDS()
        file_out.create_parameter('Charge', self.species[species_name].total_charge, 'double', parUnits='C')
        file_out.create_column('x', x, 'double', colUnits='m')
        file_out.create_column('xp', xp, 'double', colUnits='')
        file_out.create_column('y', y, 'double', colUnits='m')
        file_out.create_column('yp', yp, 'double', colUnits='')
        file_out.create_column('t', t, 'double', colUnits='s')
        file_out.create_column('p', p, 'double', colUnits='m$be$nc')
        file_out.save_sdds(file_name, dataMode='binary')
        
        return 0

    def write_opal(self, file_name, species_name):
        """Write a file to OPAL-readable format.
        :file_name: name of file to write to
        """
        coordinates = self.species[species_name].coordinates
        N = self.species[species_name].macroparticle_count
        coordinates[:, 4] -= np.mean(coordinates[:, 4])

        np.savetxt(file_name, coordinates, header='{}'.format(N), comments='')
        return 0

    def write_genesis(self, file_name, species_name, version='2.0'):
        """Write a file to genesis-readable format.
        :file_name: name of file to write to
        """
        
        # Genesis reads in external files as ASCII with the column format (per documentation):
        # X - position in x in meters
        # PX or XPRIME - momentum in x normalized to mc or divergence in x, respectively
        # Y - position in y in meters
        # PY or YPRIME - momentum in y normalized to mc or divergence in y, respectively
        # T or Z - longitudinal position in seconds or meters, respectively
        # P or GAMMA - total momentum or energy, normalized to mc or mc2, respectively.
        #
        # file has a header of the form
        # ? VERSION = 0.1
        # ? COLUMNS X PX Y PY T P
        # and the first line of data has to be the number of input particles
        
        X  = self.species[species_name].x
        PX = self.species[species_name].ux
        Y  = self.species[species_name].y
        PY = self.species[species_name].uy
        T  = self.species[species_name].ct / constants.c
        T = T - np.average(T)
        P  = self.species[species_name].pt
        
        vers_str = '? VERSION = '+version
        charge_str = '? CHARGE = '+str(abs(self.species[species_name].total_charge))
        size_str = '? SIZE = '+str(len(X))
        clmns_str = '? COLUMNS X PX Y PY T P'
        
        f = open(file_name, 'w')
        f.write(vers_str+'\n')
        f.write(charge_str+'\n')
        f.write(size_str+'\n')
        f.write(clmns_str+'\n')
        
        f.close()
        
        df = pd.DataFrame([X, PX, Y, PY, T, P]).T
        
        df.to_csv(file_name, mode='a', sep=' ', header=None, index=None)
                
        return 0

    def write(self, filename, code, species_name=None, **kwargs):
        """
        Write output file.
        :param filename: (str) Name of the file to write to.
        :param species_name: (str) Name of the species in the `species` dict to write out. Default is '{sn}'
        :param code: Name of the code format to use for writing.
        :param kwargs: Additional options to pass to the writer.
        :return:
        """.format(sn=_DEFAULT_SPECIES_NAME)

        if not species_name:
            species_name = _DEFAULT_SPECIES_NAME

        self._writers[code](filename, species_name, **kwargs)

        return filename


def _get_step_number(key, prefix='Step#'):
    result = key.split(prefix)[-1]
    return result
