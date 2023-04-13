import h5py as h5
import pandas as pd
import numpy as np
from rsbeams.rsdata.SDDS import writeSDDS, readSDDS
from rsbeams.rsptcls.species import Species
from scipy import constants

_DEFAULT_SPECIES_NAME = 'species_0'


def read_elegant(file_name: str) -> Species:
    """Read in a file from elegant output.

    elegant coordinates are as follows:
    x  -- horizontal offset, m
    x' -- horizontal angle, rad
    y  -- vertical offset, m
    y' -- vertical angle, rad
    t  -- time of flight, sec
    p  -- longitudinal momentum, mc

    Args:
        file_name: (str) Path to elegant phase space data file.

    Returns: rsbeams.rsptcls.species.Species

    """

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


def read_opal(file_name: str, step_number: int or None = None) -> Species:
    """Read in a file from OPAL output.

    opal coordinates are as follows:
    x  -- horizontal offset, m
    xp -- horizontal momentum, betax*gamma
    y  -- vertical offset, m
    yp -- vertical momentum, betay*gamma
    z  -- Position relative to ? (some sort of reference), m
    pz  -- longitudinal momentum, betaz*gamma

    Args:
        file_name: (str) Path to phase space or monitor dump.
        step_number: (int) Step number to read from phase space dump. If not given then the last step will be read.

    Returns: rsbeams.rsptcls.species.Species

    """

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
        total_charge = pcdata[loc].attrs['CHARGE']

    # TODO: This shouldn't be specific to electrons
    species = Species(particle_data, charge=-1, mass=0.511e6, total_charge=total_charge)

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

    def read(self, file_name: str, file_format: str, species_name: str or None = None, **kwargs) -> None:
        """Read particle data file.
        Readers are available for codes ([file_format specifier]: [code name]):
            elegant: elegant
            opal: OPAL
            genesis: genesis1.3

        Args:
            file_name: (str) Path to file which will be read.
            file_format: (str) Format specifier (see listing above).
            species_name: (str) [optional] Name to register species data under.
            **kwargs: (dict) Additional arguments that will be passed to the reader.

        Returns: None

        """
        reader = self._get_reader(file_format)
        if not species_name:
            species_name = _DEFAULT_SPECIES_NAME
        species = reader(file_name, **kwargs)
        self._species_insert(species_name, species)

    def write_elegant(self, file_name: str, species_name: str) -> int:
        """Write a file to elegant-readable format.

        Args:
            file_name: (str) Path for file to be written to.
            species_name: (str) Name of Species to write out.

        Returns: (int) 0 on success

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

    def write_opal(self, file_name: str, species_name: str) -> int:
        """Write a file to OPAL-readable format.

        Args:
            file_name: (str) Path for file to be written to.
            species_name: (str) Name of Species to write out.

        Returns: (int) 0 on success

        """
        coordinates = self.species[species_name].coordinates
        N = self.species[species_name].macroparticle_count

        np.savetxt(file_name, coordinates, header='{}'.format(N), comments='')
        return 0

    def write_genesis(self, file_name: str, species_name: str, version: str = '2.0') -> int:
        """Write a file to genesis-readable format.

        Args:
            file_name: (str) Path for file to be written to.
            species_name: (str) Name of Species to write out.
            version: (str) [default='2.0'] VERSION specification in header of the file.

        Returns: (int) 0 on success

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
        
        X = self.species[species_name].x
        PX = self.species[species_name].ux
        Y = self.species[species_name].y
        PY = self.species[species_name].uy
        T = self.species[species_name].ct / constants.c
        T = T - np.average(T)
        P = self.species[species_name].pt

        header = ''
        header += '? VERSION = '+version + '\n'
        header += '? CHARGE = '+str(self.species[species_name].total_charge) + '\n'
        header += '? SIZE = '+str(len(X)) + '\n'
        header += '? COLUMNS X PX Y PY T P' + '\n'

        with open(file_name, 'w') as ff:
            ff.write(header)
        
        df = pd.DataFrame([X, PX, Y, PY, T, P]).T
        
        df.to_csv(file_name, mode='a', sep=' ', header=False, index=False)
                
        return 0

    def write(self, filename: str, code: str, species_name: str or None = None, **kwargs) -> str:
        """Write output file.
        Writers are available for codes ([file_format specifier]: [code name]):
            elegant: elegant
            opal: OPAL
            genesis: genesis1.3
        
        Args:
            filename: (str) Path to write file to.
            code: (str) Simulation code to format output file for. See listing above.
            species_name: (str or None) [default=None] Name of species to write. Defaults to first registered species. 
            **kwargs: (dict) Additional arguments to be passed to the writer.

        Returns: (str) name of output file

        """.format(sn=_DEFAULT_SPECIES_NAME)

        if not species_name:
            species_name = _DEFAULT_SPECIES_NAME

        self._writers[code](filename, species_name, **kwargs)

        return filename


def _get_step_number(key, prefix='Step#'):
    result = key.split(prefix)[-1]
    return result
