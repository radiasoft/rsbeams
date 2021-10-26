from future.builtins import str
import numpy as np
from struct import pack
from sys import byteorder
from types import GeneratorType
from rsbeams.rsdata import utils
from .datatypes import supported_namelists
from .struct_data import StructData
# TODO: Would be nice to refactor the old camel case convention variables
# TODO: Add multipage write support - mostly means defining how they are input
# TODO: There may be initial nuance with the row count parameter. See no_row_counts in &data command from standard.
# TODO: Need to support additional_header_lines option (never actually seen this used though)
# TODO: Will need to test buffer read after full implementation. looks like my first tests may have led me astray.
#      multipage ascii parameter read is much slower. probably due to the need to catch comments during buffer creation
# TODO: Need to add support for SDDS versions 2-4
# TODO: Handle when readSDDS.read is called multiple times for same instance. Current behavior appends to .parameters and .columns
#       probably want to have it do nothing and print a statement to set a flag to overwrite.

# Accepted namelists in header. &data is treated as a special case.
# Types with data outside the header should be appended to the end of the list for _initialize_data_arrays
sdds_namelists = ['&associate', '&description', '&include', '&column', '&parameter', '&array']


class readSDDS:
    """
    Class for reading SDDS data files.
    Usage:
        Call `read` method to to load data from the SDDS file into memory.

    Caveats:
        - System is assumed little-endian
        - No array data (only parameters and columns)
    """

    def __init__(self, input_file, buffer=True, max_string_length=100):
        """
        Initialize the read in.

        Parameters
        ----------
        input_file: str
            Name of an SDDS file to read.
        buffer: Boolean
            If true then the file is entered into memory and closed before data is read. This may result in faster
            read times in some cases but only if the file is not on the order of available system memory.
        max_string_length: Int
            Upper bound on strings that can be read in. Should be at least as large as the biggest string in the file.
        """
        self.verbose = False
        self.buffer = buffer
        self.openf = open(input_file, 'rb')
        self.position = 0

        self.max_string_length = max_string_length
        self.header = []
        self._variable_length_records = False
        self._data_mode = None
        self.header_end_pointer = 0  #
        self._header_line_count = 1

        self._parameter_keys = []
        self.parameters = None

        self._column_keys = []
        self.columns = None

        self._array_keys = []
        self.array_size = 0
        self.arrays = None

        # Hold objects for different allowed types
        self.data = {key: [] for key in supported_namelists.keys()}

        # Read and Parse header to start
        self._read_header()
        if buffer:
            self.openf.seek(0)
            buffer = self.openf.read()
            self.openf.close()
            self.openf = buffer
        self._parse_header()
        self._initialize_data_arrays()
        self._data_descriptions = {k: utils.list_to_dict([f.fields for f in v], 'name') for k, v in self.data.items() if k != '&data'}


    @property
    def parameters(self):
        return self._parameters.data

    @parameters.setter
    def parameters(self, parameters):
        self._parameters = parameters

    @property
    def columns(self):
        return self._columns.data

    @columns.setter
    def columns(self, columns):
        self._columns = columns

    @property
    def arrays(self):
        return self._arrays.data

    @arrays.setter
    def arrays(self, arrays):
        self._arrays = arrays

    def _read_header(self):
        """
        Read in ASCII data of the header to string and organize.
        """
        if utils._read_line(self.openf).find('SDDS1') < 0:
            # First line must identify as SDDS file
            raise Exception("Header cannot be read")
        # TODO: Need a catch for &include here to at least one level of nesting
        self._header_line_count = 1
        while True:
            namelist = []
            new_line = utils._read_line(self.openf)
            self._header_line_count += 1
            if np.any([nl in new_line for nl in sdds_namelists]):
                # Log entries that describe data
                namelist.append(new_line)
                while '&end' not in new_line:
                    # Proceed reading lines until the namelist is entirely read
                    if new_line[0] == '!':
                        self._header_line_count += 1
                        continue
                    new_line = utils._read_line(self.openf)
                    namelist.append(new_line)
                    self._header_line_count += 1
                self.header.append(''.join(namelist))
            elif new_line.find('&data') == 0:
                # Final entry in the header is always &data
                namelist.append(new_line)
                while '&end' not in new_line:
                    if new_line[0] == '!':
                        self._header_line_count += 1
                        continue
                    new_line = utils._read_line(self.openf)
                    namelist.append(new_line)
                    self._header_line_count += 1
                self.header.append(''.join(namelist))
                break
            else:
                raise Exception("Header cannot be read")

        # Log data start pointer - will be adjusted to line number, if needed, after parsing
        self.header_end_pointer = self.openf.tell()

        return self.header

    def _parse_header(self):
        """
        Processes raw strings of the header into appropriate Datum structures and then creates datatype lists
        to be used by NumPy during data read.
        Returns:
            self.data holds all Datum objects
        """
        for namelist in self.header:
            namelist_type = namelist.split(maxsplit=1)[0]
            try:
                namelist_data = supported_namelists[namelist_type](namelist)
            except KeyError:
                if self.verbose:
                    print(namelist_type, namelist, 'not parsed')
                continue
            self.data[namelist_type].append(namelist_data)

        if self.data['&data'][0].fields['mode'] == 'ascii':
            self._data_mode = 'ascii'
            # need list of lines to use genfromtxt
            if self.buffer:
                self.openf = [line for line in self.openf.splitlines() if line.lstrip().find(b'!') != 0]
        else:
            self._data_mode = 'binary'
            # If binary the exact string lengths will be found dynamically and inserted here
            #self.max_string_length = '{}'

        return self.data

    def _initialize_data_arrays(self):
        for name in sdds_namelists[3:]:
            getattr(self, '_compose_'+name[1:]+'_datatypes')()
            if len(getattr(self, '_'+name[1:]+'_keys')) > 0:
                setattr(self, name[1:]+'s', StructData(getattr(self, '_'+name[1:]+'_keys'), self.max_string_length))

    def _compose_array_datatypes(self):
        pass

    def _compose_column_datatypes(self):
        """
        Creates lists of data types for all column fields
        Returns:

        """
        self._column_keys.append([])
        for col in self.data['&column']:
            if col.fields['type'] == 'string':
                if self._data_mode == 'binary':
                    self._column_keys[-1].append(('record_length', np.int32))
                    self._column_keys.append([])
                    if col.fields['field_length'] == 0:
                        self._variable_length_records = True
                    self._column_keys[-1].append((col.fields['name'],
                                                  col.type_key.format('{}')))
                else:
                    self._column_keys[-1].append((col.fields['name'],
                                                  col.type_key.format(self.max_string_length)))
            else:
                self._column_keys[-1].append((col.fields['name'], col.type_key))
        if len(self._column_keys[0]) == 0:
            # Need empty checks to succeed
            self._column_keys.pop(0)
            return
        if self.data['&data'][0].fields['mode'] == 'ascii':
            self._column_keys = [[g for f in self._column_keys for g in f]]

    def _compose_parameter_datatypes(self):
        """
        Creates a list of lists of data types for all parameter fields.
        Each sublist contains a single parameter for reading.
        Parameters are composed into a single structured array after the data read.
        Returns:

        """
        for par in self.data['&parameter']:
            if par.fields['fixed_value'] is None:
                if par.fields['type'] == 'string':
                    if self._data_mode == 'binary':
                        self._parameter_keys.append([('record_length', np.int32)])
                        self._variable_length_records = True
                        self._parameter_keys.append([(par.fields['name'],
                                                         par.type_key.format('{}'))])
                    else:
                        self._parameter_keys.append([(par.fields['name'],
                                                         par.type_key.format(self.max_string_length))])
                else:
                    self._parameter_keys.append([(par.fields['name'], par.type_key)])

        if self._data_mode == 'binary':
            # Binary always has count and it is before listed parameters start
            self._parameter_keys.insert(0, [('row_counts', utils.data_types['long'])])
        elif not self.data['&data'][0].fields['no_row_counts'] and len(self.data['&column']) > 0:
            # ASCII: count may not be included and will be at the end of the parameters
            self._parameter_keys.append([('row_counts', utils.data_types['long'])])
        else:
            if len(self._parameter_keys[0]) == 0:
                # Need empty checks to succeed
                self._parameter_keys.pop(0)

    def _get_reader(self):
        if self._data_mode == 'ascii':
            reader = np.genfromtxt
        else:
            if self.buffer:
                reader = np.frombuffer
            else:
                reader = np.fromfile

        return reader

    def _get_column_count(self, position):
        if self.data['&data'][0].fields['mode'] == 'ascii':
            count = self._get_reader()(self.openf, skip_header=position, dtype=np.int32, max_rows=1, comments='!')
        else:
            count = self._get_reader()(self.openf, offset=position, dtype=np.int32, count=1)[0]

        return count

    def _get_position(self, parameter_size, column_size, page):
        if self.buffer:
            position = (parameter_size + column_size) * page + self.header_end_pointer + self.position
        else:
            # If a file object is being passed to NumPy read tools
            #  the pointer will just always be right after the last read
            position = 0

        return position

    def _check_file_end(self, position):
        if self.buffer:
            if len(self.openf) <= position:
                return True
        else:
            # Usually will catch binary end
            pointer = self.openf.tell()

            if not self.openf.read(1):
                return True
            while True:
                line = self.openf.readline()
                if not line.strip():
                    return True
                else:
                    break

            self.openf.seek(pointer)

        return False

    def _get_ascii_row_count(self, start_position):
        line_count = 0
        if self.buffer:
            # openf is a list
            for line in self.openf[start_position:]:
                if not line:
                    return line_count
                else:
                    line_count += 1
        else:
            start_position = self.openf.tell()
            while True:
                line = self.openf.readline()
                if not line.strip():
                    self.openf.seek(start_position)
                    return line_count
                else:
                    line_count += 1

    def read(self, pages=None):
        """
        Reads all data types stored in the loaded SDDS file.
        Data is stored by field type (parameter, column, array) as attributes of readSDDS.

        Args:
            pages: If None then all pages are read. Otherwise should be an iterable object specifying
            the page numbers to be read, indexed to 0.

            e.g. pages=[0, 4, 9, 10]

        Returns:

        """
        # Always start after the header
        if not self.buffer:
            self.openf.seek(self.header_end_pointer)

        # Select pages to be stored during the read - if None then store everything
        # pages: internal counter for page numbers
        # user_pages: store what pages are returned 
        if pages:
            user_pages = pages
        else:
            user_pages = utils.iter_always()
        pages = utils.iter_always()

        # Buffer reading has no pointer, you always start at the beginning so we need to move the start offset
        if self._data_mode == 'binary':
            position = 0 + self.header_end_pointer * self.buffer
        else:
            position = 0 + self._header_line_count * self.buffer

        # TODO: Make this a function that checks if should be called
        for i, par in enumerate(self._parameter_keys):
            if par[0][0] == 'row_counts':
                rc_index = i
                break

        for page in pages:
            if not isinstance(user_pages, GeneratorType) and page > np.max(user_pages):
                print(page, user_pages)
                if self.verbose:
                    print('stopping here')
                break

            if self._check_file_end(position):
                # TODO: Get the logic right here. Should not trigger unless user requests a bad page.
                if page in user_pages:
                    if self.verbose:
                        print('Could not read page {}'.format(page))
                break

            # parameters are always read because we need to know if column_rows changes between pages
            parameter_data, position = self._get_parameter_data(self._parameter_keys, position)
            # save data if needed
            if parameter_data and (page in user_pages):
                self._parameters.add(parameter_data)

            if len(self.data['&column']) == 0:
                row_count = 0
            elif not self.data['&data'][0].fields['no_row_counts']:
                row_count = utils.get_entry_from_parameters(parameter_data, 'row_counts')  #parameter_data[rc_index][0][0][0]
            else:
                row_count = self._get_ascii_row_count(position)

            if row_count is 0:
                continue
            if (isinstance(user_pages, GeneratorType) or (page in user_pages)) or self._variable_length_records:
                column_data, position = self._get_column_data(self._column_keys, position, row_count)
                self._columns.add(column_data)
            else:
                # still need to update position what would have been read
                position += np.dtype(self._column_keys[0]).itemsize * row_count
        if self._parameters:
            self._parameters.concat()
        if self._columns:
            self._columns.concat()

    def _get_parameter_data(self, data_keys, position):
        data_arrays = [[]]
        for dk in data_keys:
            if self._variable_length_records:
                try:
                    record_length = data_arrays[0][-1]['record_length']
                    dk = [(dk[0][0], dk[0][1].format(record_length[0]))]
                except (ValueError, IndexError):
                    pass
            if self.data['&data'][0].fields['mode'] == 'ascii':
                new_array = self._get_reader()(self.openf, skip_header=position, dtype=dk, max_rows=1,
                                               comments='!', deletechars='', unpack=True, delimiter='\n\t')
                if self.buffer:
                    position += 1
            else:
                new_array = self._get_reader()(self.openf, dtype=dk, count=1, offset=position)
                if self.buffer:
                    position += np.dtype(dk).itemsize
            data_arrays[-1].append(new_array)

        return data_arrays, position

    def _get_column_data(self, data_keys, position, row_count):
        # TODO: Account for now_row_count possibility
        data_arrays = []
        # Hand variable record length read by iterating through rows
        if len(data_keys) > 1:
            for row in range(row_count):
                data_arrays.append([])
                for i, dk in enumerate(data_keys):
                    if self._variable_length_records:
                        dk = dk.copy()
                        try:
                            record_length = data_arrays[-1][-1]['record_length']
                            dk[0] = (dk[0][0], dk[0][1].format(record_length[0]))
                        except (ValueError, IndexError):
                            pass
                    if self.data['&data'][0].fields['mode'] == 'ascii':
                        new_array = self._get_reader()(self.openf, skip_header=position, dtype=dk, max_rows=1,
                                                       comments='!', deletechars='')
                        if self.buffer:
                            position += 1
                    else:
                        new_array = self._get_reader()(self.openf, dtype=dk, count=1, offset=position)
                        if self.buffer:
                            position += np.dtype(dk).itemsize
                    data_arrays[-1].append(new_array)
        else:
            # If no variable records all rows can be read at once
            data_arrays.append([])
            dk = data_keys[0]
            if self.data['&data'][0].fields['mode'] == 'ascii':
                new_array = self._get_reader()(self.openf, skip_header=position, dtype=dk, max_rows=row_count,
                                               comments='!', deletechars='')
                if self.buffer:
                    position += 1 * ((0 if not row_count else row_count) + 1 * int(self.data['&data'][0].fields['no_row_counts']))
            else:
                new_array = self._get_reader()(self.openf, dtype=dk, count=row_count, offset=position)
                if self.buffer:
                    position += np.dtype(dk).itemsize * row_count
            data_arrays[-1].append(new_array)

        return data_arrays, position

    def parameter_symbol(self, parameter):
        return self._data_descriptions['&parameter'][parameter]['symbol']

    def parameter_units(self, parameter):
        return self._data_descriptions['&parameter'][parameter]['units']

    def parameter_description(self, parameter):
        return self._data_descriptions['&parameter'][parameter]['description']

    def parameter_format_string(self, parameter):
        return self._data_descriptions['&parameter'][parameter]['format_string']

    def parameter_type(self, parameter):
        return self._data_descriptions['&parameter'][parameter]['type']

    def parameter_fixed_value(self, parameter):
        return self._data_descriptions['&parameter'][parameter]['fixed_value']

    def column_symbol(self, column):
        return self._data_descriptions['&column'][column]['symbol']

    def column_units(self, column):
        return self._data_descriptions['&column'][column]['units']

    def column_description(self, column):
        return self._data_descriptions['&column'][column]['description']

    def column_format_string(self, column):
        return self._data_descriptions['&column'][column]['format_string']

    def column_type(self, column):
        return self._data_descriptions['&column'][column]['type']

headSDDS = "SDDS1\n"
columnAttributeStr = {'colName': ['name=', '{}'], 'colType': ['type=', '{}'], 'colUnits': ['units=', '"{}"'],
                      'colSymbol': ['symbol=', '"{}"'], 'colFormatStr': ['format_string=', '"{}"'],
                      'colDescription': ['description=', '"{}"'], 'colFieldLen': ['field_length=', '"{}"']}
parameterAttributeStr = {'parName': ['name=', '{}'], 'parType': ['type=', '{}'], 'parUnits': ['units=', '"{}"'],
                         'parSymbol': ['symbol=', '"{}"'], 'parFormatStr': ['format_string=', '"{}"'],
                         'parDescription': ['description=', '"{}"'], 'parFixedVal': ['fixed_value=', '{}']}


class writeSDDS:
    """
    Implements an SDDS class in Python.
    Can be used to write out data stored as NumPy arrays or single values stored as a variable
    Included methods:
        SDDS.create_column
        SDDS.create_param
        SDDS.save_sdds
    Does not support creating multi-page SDDS files at this time.
    Acceptable values for colType/parType:
        short
        long
        double
        character (not recommended)
        string    (not recommended)

    @author: Chris
    """

    # TODO: Test binary writeout more carefully. (Multiple types)

    sddsIdentifier = headSDDS

    key_indentity = {'double': 'd', 'short': 's', 'long': 'i', 'string': None}

    def __init__(self, page=1, readInFormat='numpyArray'):
        """
        Initialize SDDS object for storing parameter/column data and writing.
        """
        self.format = readInFormat
        self.endianness = byteorder
        self.page = page
        self.columns = []
        self.parameters = []
        self.column_key = '='
        self.dataMode = 'ascii'

    def create_column(self, colName, colData, colType,
                      colUnits='', colSymbol='', colFormatStr='', colDescription='', colFieldLen=0):
        """
        Creates a column data object that can be written to file.

        Parameters
        ----------
        colName: str
            Name of the column.
        colData: ndarray (Data type must match 'colType')
        colType: str
            Data type for the column. Must match data type contained in 'colData'. See Description of the class
            for available data types to write.
        colUnits: str (optional)
            String with the units of the column. To be written out to file.
        colSymbol: str (optional)
            Optional symbol string that can be written out to the file. See SDDS manual for syntax.
        colFormatStr: str (optional)
            May specify the form of the printf string for use by SDDS.
        colDescription: str (optional)
            Optional description of the column to write to file.
        colFieldLen: int (optional)
            For ASCII data, the optional field length field specifies the number of characters occupied by
            the data for the column. If zero, the data is assumed to be bounded by whitespace characters. If
            negative, the absolute value is taken as the field length, but leading and trailing whitespace characters
            will be deleted from string data. This feature permits reading fixed-field-length FORTRAN
            output without modification of the data to include separators
        """
        # TODO: Need to implement handler for FieldLen parameter setting once I support strings
        assert colType in self.key_indentity.keys(), "Allowed data types are: 'double', 'short', 'long' "

        column_data = {'colName': colName, 'colData': colData, 'colType': colType, 'colUnits': colUnits,
                       'colSymbol': colSymbol, 'colFormatStr': colFormatStr, 'colDescription': colDescription,
                       'colFieldLen': colFieldLen}

        self.columns.append(column_data)

    def create_parameter(self, parName, parData, parType,
                         parUnits='', parSymbol='', parFormatStr='', parDescription='', parFixedVal=None):
        """
        Creates a parameter data object that can be written to file.

        Parameters
        ----------
        parName: str
            Name of the parameter.
        parData: short, long, or double
            Data being written to the SDDS file.
        parType: str
            Data type for the parameter. Must match data type for the variable being written.
            See Description of the class for available data types to write.
        parUnits: str (optional)
            String with the units of the parameter. To be written to file.
        parSymbol: str (optional)
            Optional symbol string that can be written out to the file. See SDDS manual for syntax.
        parFormatStr: str (optional)
            May specify form of the printf string for use by SDDS.
        parDescription: str (optional)
            Optional description of the parameter to be written to the file.
        parFixedVal: (float, int, str) (optional) May be used to store fixed data in the file header, instead
            of the data body. If given then parData must be None.
        """
        assert parType in self.key_indentity.keys(), "Allowed data types are: 'double', 'short', 'long' "

        # Handle variable data type possibility for fixed_value
        if parFixedVal:
            assert parData is None, "Using fixed value requires parData=None"
            if type(parFixedVal) in [float, int]:
                pass
            elif type(parFixedVal) == str:
                parFixedVal = '"{}"'.format(parFixedVal)

        parameter_data = {'parName': parName, 'parData': parData, 'parType': parType, 'parUnits': parUnits,
                          'parSymbol': parSymbol, 'parFormatStr': parFormatStr, 'parDescription': parDescription,
                          'parFixedVal': parFixedVal}

        self.parameters.append(parameter_data)

    def _write_header(self, outputFile):
        columnString = ''
        parameterString = ''
        outputFile.write(self.sddsIdentifier.encode())
        outputFile.write('!# {}-endian\n'.format(self.endianness).encode())

        for parameter in self.parameters:
            for key, value in parameter.items():
                if (key in parameterAttributeStr.keys()) and value:
                    parameterString = ''.join((parameterString, '{field}{value}, '.format(field=parameterAttributeStr[key][0],
                                                                                          value=parameterAttributeStr[key][1].format(value))))
            outputFile.write('&parameter {}&end\n'.format(parameterString).encode())
            parameterString = ''

        for column in self.columns:
            for key, value in column.items():
                if key in columnAttributeStr:
                    # split logic to avoid the hassle of testing existance of arrays separately
                    if value:
                        columnString = ''.join((columnString, '{field}{value}, '.format(field=columnAttributeStr[key][0],
                                                                    value=columnAttributeStr[key][1].format(value))))
            outputFile.write('&column {} &end\n'.format(columnString).encode())
            columnString = ''

        outputFile.write('&data mode={}, &end\n'.format(self.dataMode).encode())

    def save_sdds(self, fileName, dataMode='ascii'):
        """
        Saves the parameters and columns to file. Parameters and columns are written to the file in the order
        that they were created in the writeSDDS object.

        Parameters
        ----------
        fileName: str
            Name of the file to be written.
        dataMode: Either 'ascii' or 'binary'
            Write mode for the file. Data will either be written to the file in ascii or binary mode.
        """

        self.dataMode = dataMode
        outputFile = open(fileName, 'wb')

        # Verify Column Data Integrity
        if len(self.columns) > 1:
            try:
                column_data = np.column_stack([columns['colData'] for columns in self.columns])
            except ValueError:
                print('ERROR: All columns on a page must have same length')
                return
        elif len(self.columns) == 1:
            column_data = self.columns[0]['colData']
        else:
            column_data = np.empty([0])

        # Begin header writeout
        self._write_header(outputFile)

        # Write row count. Write 0 if no rows.
        if self.dataMode == 'binary':
            # Row count precedes parameter entries in a binary file
            outputFile.write(pack('I', column_data.shape[0]))

        # Write Parameters
        for parameter in self.parameters:
            # Pass if fixed_value used (parData will be None)
            if parameter['parData'] is not None:
                if self.dataMode == 'ascii':
                        outputFile.write('{}\n'.format(parameter['parData']).encode())
                if self.dataMode == 'binary':
                        outputFile.write(pack('={}'.format(self.key_indentity[parameter['parType']]),
                                              parameter['parData']))

        # Write row count. Write 0 if no rows.
        if self.dataMode == 'ascii':
            # Row count follows parameter entries in an ascii file
            # Should not appear in ascii if 0
            if column_data.shape[0]:
                outputFile.write('{}\n'.format(column_data.shape[0]).encode())

        # Write Columns
        if self.dataMode == 'ascii':
            np.savetxt(outputFile, column_data)
        elif self.dataMode == 'binary':
            column_data.tofile(outputFile)
            # outputFile.write(pack('=i7si8s', *column_data))  Fixed Test for ascii output Chicago New York
        else:
            print("NOT A DEFINED DATA TYPE")

        outputFile.close()
