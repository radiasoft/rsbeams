from future.builtins import str
import numpy as np
import shlex
from struct import pack, unpack, calcsize, error
from sys import byteorder
from copy import copy
# TODO: Will probably need a new dict object that uses ordered Dict for py2 and userdict or standard dict for py3
# TODO: Would be nice to refactor the old camel case convention variables
# TODO: Add multipage write support - mostly means defining how they are input
# TODO: There may be initial nuance with the row count parameter. See no_row_counts in &data command from standard.
# TODO: Need to support additional_header_lines option (never actually seen this used though)


def _return_type(string):
    target = 'type='
    type_field_position = string.find(target)
    if type_field_position < 0:
        return None
    for key in data_types.keys():
        if string[type_field_position:type_field_position+len(key)].find(key) > -1:
            return key
    return None


def _read_line(open_file):
    line = ''
    while line == '\n' or line == '' or line[0] == '!':
        line = str(open_file.readline(), 'latin-1')
    return line


def _shlex_split(line):
    # Cannot just use shlex.split
    # Need to ignore single-quotes due to their use as common characters in accelerator notation
    lex = shlex.shlex(line, posix=True)
    lex.quotes = '"'
    lex.whitespace_split = True
    lex.commenters = ''

    return list(lex)


def iter_always():
    i = 0
    while True:
        yield i
        i += 1


# Accepted namelists in header. &data is treated as a special case.
sdds_namelists = ['&column', '&parameter', '&description', '&array', '&include']
data_types = {'double': np.float64, 'short': np.int32, 'long': np.int64, 'string': 'S{}', 'char': np.char}


class Datum:
    def __init__(self, namelist):
        self.namelist = namelist
        self.type_key = None
        self.parse_namelist()

    def parse_namelist(self):
        namelist = _shlex_split(self.namelist)
        for name in self.fields.keys():
            for entry in namelist:
                if name in entry[:len(name)]:
                    start = entry.find('=')
                    self.fields[name] = entry[start + 1:].rstrip(',')

    def set_data_type(self):
        data_type = self.fields['type']
        self.type_key = data_types[data_type]


class Parameter(Datum):
    def __init__(self, namelist):
        self.fields = {'name': '', 'symbol': '', 'units': '', 'description': '',
                       'format_string': '', 'type': '', 'fixed_value': None}
        super().__init__(namelist)

        if self.fields['fixed_value'] is not None:
            self.parse_fixed_value()
        else:
            self.set_data_type()

    def parse_fixed_value(self):
        fixed_value = self.fields['fixed_value']
        data_type = self.fields['type']
        if fixed_value:
            if data_type in ['short', 'long']:
                self.fields['fixed_value'] = int(fixed_value)
            elif data_type == 'double':
                self.fields['fixed_value'] = float(fixed_value)
            elif data_type in ['string', 'char']:
                try:
                    self.fields['fixed_value'] = str(fixed_value).decode('utf-8')  # TODO: Look into this warning
                except AttributeError:
                    self.fields['fixed_value'] = str(fixed_value)


class Column(Datum):
    def __init__(self, namelist):
        self.fields = {'name': '', 'symbol': '', 'units': '', 'description': '',
                       'format_string': '', 'type': '', 'field_length': 0}
        super().__init__(namelist)
        self.set_data_type()

    def set_data_type(self):
        # Override to account for field_length setting
        data_type = self.fields['type']
        if data_type == 'string':
            if self.fields['field_length'] == 0:
                self.type_key = data_types[data_type]
            else:
                self.type_key = data_types[data_type].format(np.abs(self.fields['field_length']))
        else:
            self.type_key = data_types[data_type]


class Data(Datum):
    def __init__(self, namelist):
        self.fields = {'mode': 'binary', 'lines_per_row': 1, 'no_row_counts': 0, 'additional_header_lines': 0}
        super().__init__(namelist)


# Available namelists from SDDS standard
supported_namelists = {'&parameter': Parameter, '&column': Column, '&data': Data}


class readSDDS:
    """
    Class for reading SDDS data files.
    Usage:
        Call `read` method to to load data from the SDDS file into memory.

    Caveats:
        - System is assumed little-endian
        - Data stored in binary
        - No array data (only parameters and columns)
        - Strings stored as fixed are not included in parameters (strings stored in binary format are included)
        - Files that store string data in columns are not currently supported
    """

    def __init__(self, input_file, verbose=False, buffer=True, max_string_length=100):
        """
        Initialize the read in.

        Parameters
        ----------
        input_file: str
            Name of binary SDDS file to read.
        verbose: Boolean
            Print additional data about detailing intermediate read in process.
        max_string_length: Int
            Upper bound on strings that can be read in. Only used for formatting read from ASCII files.
        """
        self._input_file = input_file
        self.verbose = verbose
        self.buffer = buffer
        if buffer:
            openf = open(input_file, 'rb')
            self.openf = openf.read()
            openf.close()
        else:
            self.openf = open(input_file, 'rb')

        self.max_string_length = max_string_length
        self.header = []
        self._variable_length_records = False

        self.param_key = ['=i']  # Include row count with parameters
        self.column_key = '='
        self.header_end_pointer = 0  #

        self._parameter_keys = []
        self.param_size = 0
        self.param_names = ['rowCount']
        self.parameters = {}

        self._column_keys = []
        self.column_size = 0
        self.column_names = {}
        self.columns = False

        # Hold objects for different allowed types
        self.data = {key: [] for key in supported_namelists.keys()}

        # Read and Parse header to start
        self._read_header()
        self._parse_header()

    def _read_header(self):
        """
        Read in ASCII data of the header to string and organize.
        """
        if _read_line(self.openf).find('SDDS1') < 0:
            # First line must identify as SDDS file
            raise Exception("Header cannot be read")
        # TODO: Need a catch for &include here to at least one level of nesting
        while True:
            namelist = []
            new_line = _read_line(self.openf)
            if np.any([nl in new_line for nl in sdds_namelists]):
                # Log entries that describe data
                namelist.append(new_line)
                while '&end' not in new_line:
                    # Proceed reading lines until the namelist is entirely read
                    if new_line[0] == '!':
                        continue
                    new_line = _read_line(self.openf)
                    namelist.append(new_line)
                self.header.append(''.join(namelist))
            elif new_line.find('&data') == 0:
                # Final entry in the header is always &data
                namelist.append(new_line)
                while '&end' not in new_line:
                    if new_line[0] == '!':
                        continue
                    new_line = _read_line(self.openf)
                self.header.append(''.join(namelist))
                break
            else:
                raise Exception("Header cannot be read")

        # Log data start position
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
                print(namelist_type, namelist, 'not parsed')  # TODO: TEMP printout
                continue
            self.data[namelist_type].append(namelist_data)

        if self.data['&data'][0].fields['mode'] == 'binary':  # TODO: Is there any scenario where this more than one &data field?
            self.max_string_length = '{}'

        return self.data

    def _compose_column_datatypes(self):
        """
        Creates lists of data types for all column fields
        Returns:

        """
        self._column_keys.append([])
        for col in self.data['&column']:
            if col.fields['type'] == 'string':
                self._column_keys[-1].append(('record_length', np.int32))
                self._column_keys.append([])
                self._column_keys[-1].append((col.fields['name'],
                                              col.type_key.format(self.max_string_length)))

                if col.fields['field_length'] == 0:
                    self._variable_length_records = True
            else:
                self._column_keys[-1].append((col.fields['name'], col.type_key))

        if self.data['&data'][0].fields['mode'] == 'ascii':
            self._column_keys = [g for f in self._column_keys for g in f]

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
                    self._parameter_keys.append([('record_length', np.int32)])
                    self._parameter_keys.append([(par.fields['name'],
                                                     par.type_key.format(self.max_string_length))])
                    self._variable_length_records = True
                else:
                    self._parameter_keys.append([(par.fields['name'], par.type_key)])

    def _get_reader(self):
        if self.data['&data'].fields['mode'] == 'ascii':
            reader = np.genfromtxt
        else:
            if self.buffer:
                reader = np.frombuffer
            else:
                reader = np.fromfile

        return reader

    def _get_column_count(self, position):
        if self.data['&data'].fields['mode'] == 'ascii':
            count = self._get_reader()(self.openf, skip_header=position, dtype=np.int32, max_rows=1, comments='!')
        else:
            count = self._get_reader()(self.openf, offset=position, dtype=np.int32, count=1)[0]

        return count

    def _get_position(self, parameter_size, column_size, page):
        # TODO: Pages indexed to 0
        if self.buffer:
            position = (parameter_size + column_size) * page + self.header_end_pointer + self.position
        else:
            # If a file object is being passed to NumPy read tools
            #  the pointer will just always be right after the last read
            position = 0

        return position

    def _check_file_end(self, position):
        if self.buffer:
            if len(self.openf) == position:
                return True
        else:
            pointer = self.openf.tell()
            if not self.openf.read(1):
                return True
            self.openf.seek(pointer)

        return False

    def read2(self, pages):
        if pages:
            user_pages = pages
        else:
            user_pages = iter_always()
        pages = iter_always()
        parameter_size, column_size = 0, 0  # Total size consumed so far from the file
        position = self.header_end_pointer
        # if len(self._column_keys) == 1:
        #     row_size = np.dtype(self._column_keys[0]).itemsize
        # else:
        #     row_size = -1
        #     # If there are variable records we must read every page, only pages user requests will be stored though
        #     pages = iter_always()
        # if len(self._parameter_keys) == 1:
        #     parameter_size = np.dtype(self._parameter_keys[0]).itemsize
        # else:
        #     parameter_size = -1
        #     # If there are variable records we must read every page, only pages user requests will be stored though
        #     pages = iter_always()

        for page in pages:
            # get position of the page start. if page = 1 then don't need par and col sizes anyway
            # if page > 1 then we will have the last calculated sizes. get_position will need to store the accumulated
            # value though. could be iterable.
            position = self._get_position(parameter_size, column_size, page)
            if not self.data['&data'].fields['no_row_count']:
                row_count = self._get_column_count(position)
            # if self._check_file_end(position):
            #     # TODO: May not need this if column_count catches file end
            #     if page in user_pages:
            #         print('Could not read page {}'.format(page))
            #     break

            # based on position and data key get the data and update the parameter_size if it was unknown (<0)
            parameter_data, parameter_size = _get_parameter_data(self._parameter_keys, position)
            # save data if needed
            if page in user_pages:
                self.parameter_data.add(parameter_data)
            # same but for columns
            column_data, column_size = _get_column_data(self._column_keys, position + parameter_size, row_count)
            if page in user_pages:
                self.column_data.add(column_data)
            self.position = position

    def _get_parameter_data(self, data_keys, position):
        data_arrays = []
        if len(data_keys) > 1:
            for dk in data_keys:
                try:
                    record_length = data_arrays[-1]['record_length']
                    dk[0] = (dk[0][0], dk[0][1].format(record_length[0]))
                except (ValueError, IndexError):
                    pass
                if self.data['&data'].fields['mode'] == 'ascii':
                    new_array = self._get_reader()(self.openf, skip_header=position, dtype=dk, max_rows=len(dk),
                                                   comments='!')
                else:
                    new_array = self._get_reader()(self.openf, dtype=dk, count=1, offset=position)
                # print(new_array, new_array.dtype)
                data_arrays.append(new_array)
                # TODO: Would it be faster to consume the buffer as we go?
                position += np.dtype(dt).itemsize



    def _read_params(self):
        """
        Read parameter data from the SDDS file.

        Returns
        -------
        parameters: dictionary
            Dictionary object with parameters names and values.
        """

        param_data = ()
        for key in self.data['&Parameter']:
            if key[0] == 'z':
                str_length = unpack('i', self.openf.read(4))[0]
                str_size = calcsize('=' + 'c' * str_length)
                value = unpack('=' + 'c' * str_length, self.openf.read(str_size))
                value = ''.join(value)
                param_data = param_data + (value,)
            else:
                value = unpack(key, self.openf.read(calcsize(key)))
                if self.verbose:
                    print(value)
                param_data = param_data + value

        for param, value in zip(self.param_names, param_data):
            self.parameters[param] = value

        self.row_count = self.parameters['rowCount']

        return self.parameters

    def _read_columns(self):
        """
        Read column data from the SDDS file.

        Returns
        -------
        columns: ndarray
            NumPy array with column data.
        """
        # if self.columns:
        #     return np.asarray(self.columns)
        # else:
        #     pass

        try:
            self.row_count
        except AttributeError:
            self._read_params()

        self.columns = []
        # TODO: Can probably use numpy to read in columns much more quickly and split after
        for i in range(self.row_count):
            try:
                self.columns.append(unpack(self.column_key, self.openf.read(self.column_size)))
            except error:
                break

        return np.asarray(self.columns)

    def _get_row_count(self, here=False):
        """
        Get row count on a page. Will leave file position position at start of the row count that was read.
        Args:
            here (boolean): If False then get the row count of the first page. Else try to read at current position.

        Returns:

        """

        if not here:
            self.openf.seek(self.header_end_pointer)
        row_count = unpack('i', self.openf.read(calcsize('i')))[0]

        # Set file position appropriately
        if not here:
            self.openf.seek(self.header_end_pointer)
        else:
            self.openf.seek(self.openf.tell() - calcsize('i'))

        return row_count

    def _set_all_row_counts(self):
        row_count = [self._get_row_count(here=False)]
        while True:
            # Move to next page
            current_position = self.openf.tell()
            page_size = self.param_size + row_count[-1] * self.column_size
            self.openf.seek(current_position + page_size)

            # Catch the end of file
            try:
                row_count.append(self._get_row_count(here=True))
            except error as e:
                break

        self._row_count = row_count

    def read(self, pages=None):
        """
        Read page(s) from the SDDS file into memory.
        Args:
            pages (int or list or tuple): If int then the number of pages to read. If iterable then a list of
            the pages to read using 0-based indexing. If None then all pages will be read.

        Returns:

        """
        parameters = []
        columns = []

        # Construct list of pages
        if pages:
            if type(pages) == int:
                pages = np.arange(pages)
            elif type(pages) == list or type(pages) == tuple:
                pages = np.array(pages)
            else:
                raise TypeError("pages must be an int, tuple, or list")
        else:
            def iter_always():
                i = 0
                while True:
                    yield i
                    i += 1

            pages = iter_always()

        # Get sizes to skip pages if needed
        self._set_all_row_counts()
        header_size = self.header_end_pointer
        param_size = self.param_size
        # Need the size of all rows and columns (not just one row) for this
        col_sizes = self.column_size * np.array(self._row_count)

        # Read all requested pages
        for page in pages:

            # check if pages were skipped and update position
            expected_position = page * param_size + np.sum(col_sizes[:page]) + header_size
            if expected_position != self.openf.tell():
                self.openf.seek(expected_position)

            # Catch the end of file
            try:
                self._get_row_count(here=True)
            except error as e:
                # Catch if user requested an unreadable page
                if type(pages) == np.ndarray:
                    print("WARNING: Could not read page {}".format(page))
                break

            params = self._read_params()
            parameters.append(copy(params))

            if self.row_count > 0:
                cols = self._read_columns()
                columns.append(cols)

        self.parameters = parameters
        self.columns = np.asarray(columns)

        return parameters, np.asarray(columns)

    def close(self):
        # TODO: If we keep this it means people might want to call read multiple times and then pointer position may need to be reset
        """
        Close opened SDDS file.
        Returns:

        """
        self.openf.close()


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
            if parameter['parData']:
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


if __name__ == '__main__':
    ff = readSDDS('run_setup.output.sdds')
    print(ff.param_key[0], ff.column_key)
    p, c = ff.read(pages=[0, 1, 2, 4])
    print(p, c.shape)