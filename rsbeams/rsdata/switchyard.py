# Import the relevant data formats
# TODO: readers and writers should just be functions
# TODO: elegant reader can use readSDDS now
# TODO: wrap species handling after it is instantiated
# TODO: create centralized record of unit conversions by code

import pandas as pd
import numpy as np
import h5py as h5
import os, re
from scipy import constants
from rsbeams.rsptcls.species import Species
from subprocess import Popen, PIPE
from rsbeams.rsdata.SDDS import writeSDDS
from openpmd_viewer import OpenPMDTimeSeries

supported_codes = ['genesis', 'elegant', 'opal', 'warp']


class Switchyard:
    """Class for writing a particle species data from a code output to a universal format, or vice verse. 
    Can be used to load output from one code into another, or for universal data visualization."""
    
    def __init__(self, input_file, input_format):
        
        # conventions
        # x  -- horizontal position, m
        # ux -- horizontal velocity, beta_x gamma
        # y  -- vertical position, m
        # uy -- vertical velocity, beta_y gamma
        # ct -- time of flight, m
        # pt -- total momentum, beta gamma
        self.input_file = input_file
        self.supported_codes = supported_codes
        self.species = {}
        assert input_format in self.supported_codes, "{} is not supported".format(input_format)
        self.input_format = input_format
        self._readers = {'elegant': self.read_elegant, 'opal': self.read_opal, 'warp': self.read_openpmd}
        self._writers = {'elegant': self.write_elegant, 'genesis': self.write_genesis}

        self._get_reader()(file_name=input_file)
    
    def is_supported(self, code_name):
        
        if code_name in self.supported_codes:
            return True
        else:
            return False

    def supported_codes(self):
        
        return self.supported_codes

    def _get_reader(self, file_format=None):
        if not file_format:
            file_format = self.input_format
            return self._readers[file_format]

    def read_elegant(self, file_name, species_name = 'Species'):
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

        read_particle_data = 'sdds2stream -col=x,xp,y,yp,t,p {file}'.format(file=file_name)
        read_charge = 'sdds2stream -par=Charge {file}'.format(file=file_name)
        get_particle_data = Popen(read_particle_data, stdout=PIPE, stderr=PIPE, shell=True)
        get_charge_data = Popen(read_charge, stdout=PIPE, stderr=PIPE, shell=True)
        particle_data, err = get_particle_data.communicate()
        if err:
            return err
        else:
            particle_data = np.fromstring(particle_data, dtype=float, count=-1, sep=' \n').reshape(-1, 6)
        charge_data, err = get_charge_data.communicate()
        if err:
            return err
        else:
            charge_data = np.fromstring(charge_data, dtype=float, count=1, sep=' \n')[0]
            
        if species_name == 'Species':
            spec_name = species_name+'_'+str(len(self.species.keys()))
        else:
            spec_name = species_name
        self.species[spec_name] = Species(particle_data, charge=-1, mass=0.511e6, total_charge=charge_data)
        self.species[spec_name].convert_from_elegant()

        return 0
    
    def read_opal(self, file_name, step_number=None, species_name = 'Species'):
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
        # TODO: We really need to be able to pull from screens too but we'll settle for standard distribution output for now
        # TODO: We don't currently handle particle specific weight/charge in Species
        
        with h5.File(file_name, 'r') as pcdata:
            if not step_number:
                step_number = len(pcdata) - 1
            loc = 'Step#{}'.format(step_number)
            mp_count = pcdata[loc+'/z'].shape[0]
            particle_data = np.empty((mp_count, 6))
            for i, coord in enumerate(['x', 'px', 'y', 'py', 'z', 'pz']):
                print(loc+'/'+coord)
                particle_data[:, i] = pcdata[loc+'/'+coord]
            total_charge = pcdata[loc].attrs['CHARGE']
        
        # TODO: This needs to be a function --- it is used many times
        #  probably just make a species instantiation function
        if species_name == 'Species':
            spec_name = species_name+'_'+str(len(self.species.keys()))
        else:
            spec_name = species_name

        # TODO: This shouldn't be specific to electrons
        self.species[spec_name] = Species(particle_data, charge=-1, mass=0.511e6, total_charge=total_charge)
        
        return 0

    def read_genesis(self, file_name):
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

    def read_openpmd(self, file_name, species_name):
        # OpenPMDTimeSeries needs a directory and will open all OpenPMD files in that directory
        directory, fname = os.path.split(file_name)
        requested_iteration = int(re.search(r"\d+", fname).group())

        ts = OpenPMDTimeSeries(directory)
        assert requested_iteration in ts.iterations, f"Could not find iteration {requested_iteration}"
        particle_data = ts.get_particle(var_list=['x', 'y', 'z', 'ux', 'uy', 'uz', 'charge', 'mass', 'w'],
                                        pecies=species_name, iteration=requested_iteration)

        total_charge = np.sum(particle_data[:, 6] * particle_data[:, 8])
        self.species[species_name] = Species(particle_data[:, :6],
                                          charge=particle_data[:, 6],
                                          mass=particle_data[:, 7], total_charge=total_charge)
        self.species[species_name].convert_from_openpmd()
    
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

    def write_opal(self, file_name):
        """Write a file to OPAL-readable format.
        :file_name: name of file to write to
        """
        
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
        charge_str = '? CHARGE = '+str(self.species[species_name].total_charge)
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

    def write(self, filename, code, species_name='Species_0', **kwargs):
        """
        Write output file.
        :param filename: (str) Name of the file to write to.
        :param species_name: (str) Name of the species in the `species` dict to write out. Default is 'Species_0'
        :param code: Name of the code format to use for writing.
        :param kwargs: Additional options to pass to the writer.
        :return:
        """

        self._writers[code](filename, species_name, **kwargs)

        return filename