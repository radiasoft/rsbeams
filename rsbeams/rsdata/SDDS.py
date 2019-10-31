from future.builtins import str
import numpy as np
from struct import pack, unpack, calcsize, error
from sys import byteorder
from copy import copy
# TODO: Will probably need a new dict object that uses ordered Dict for py2 and userdict or standard dict for py3
# TODO: Would be nice to refactor the old camel case convention variables
# TODO: Add multipage write support - mostly means defining how they are input
# TODO: There may be initial nuance with the row count parameter. See no_row_counts in &data command from standard.


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

    def __init__(self, input_file, verbose=False):
        """
        Initialize the read in.

        Parameters
        ----------
        input_file: str
            Name of binary SDDS file to read.
        verbose: Boolean
            Print additional data about detailing intermediate read in process.
        """

        self.openf = open(input_file, 'rb')
        self.verbose = verbose
        self.header = []

        self.param_key = ['=i']  # Include row count with parameters
        self.column_key = '='
        self.header_end_pointer = 0  #
        self.string_in_params = []
        self.string_in_columns = []

        self.param_size = 0
        self.column_size = 0

        self.param_names = ['rowCount']
        self.column_names = {}
        self.parameters = False
        self.columns = False

        # Read and Parse header to start
        self._read_header()
        self._parse_header()

    def _read_header(self):
        """
        Read in ASCII data of the header to string and organize.
        """

        while True:
            new_line = str(self.openf.readline(), 'latin-1')

            if new_line.find('&data') == 0:
                self.header.append(new_line)
                break
            else:
                self.header.append(new_line)

        self.header_end_pointer = self.openf.tell()

        return self.header

    def _parse_header(self):
        # TODO: figure out why param_key is list with just the key and col_key is just the string defining to key
        """
        Parse header data to instruct unpacking procedure.
        """

        self.parsef = True

        params = []
        columns = []
        parameter_position = 0

        # Find Parameters and Column
        for line in self.header:
            if line.find('&parameter') == 0:
                params.append(line)
            if line.find('&column') == 0:
                columns.append(line)

        # Construct format string for parameters and columns
        for param in params:
            if param.find('type=string') > -1:
                if param.find('fixed_value') > -1:  # Fixed value not present means string is in binary data
                    if self.verbose:
                        print('passed')
                    pass
                else:
                    if self.verbose:
                        print('used')
                    self.param_key.append('zi')
                    parameter_position += 2
                    self.param_key.append('=')
            elif param.find('type=double') > -1:
                self.param_key[parameter_position] += 'd'
            elif param.find('type=long') > -1:
                self.param_key[parameter_position] += 'i'
            elif param.find('type=short') > -1:
                self.param_key[parameter_position] += 's'
            else:
                pass

        if self.param_key[-1] == '=':
            self.param_key.pop(-1)  # Remove the last '=' that will be added if final entry is string

        for column in columns:
            if column.find('type=double') > -1:
                self.column_key += 'd'
            elif column.find('type=long') > -1:
                self.column_key += 'i'
            elif column.find('type=short') > -1:
                self.column_key += 's'
            else:
                pass

        for param in params:
            if param.find('type=string') > -1 and param.find('fixed_value=') > -1:
                pass
            else:
                i0 = param.find('name') + 5
                ie = param[i0:].find(',')
                self.param_names.append(param[i0:i0 + ie])
        for i, column in enumerate(columns):
            i0 = column.find('name') + 5
            ie = column[i0:].find(',')
            self.column_names[column[i0:i0 + ie]] = i

        self.param_size = calcsize(self.param_key[0])
        self.column_size = calcsize(self.column_key)

        if self.verbose:
            print("Parameter unpack size: %s bytes \nColumn unpack size: %s bytes".format(
                (self.param_size, self.column_size)))

        if self.verbose:
            print("Parameter key string: %s \nColumn key string: %s".format(self.param_key, self.column_key))

        return self.param_key, self.column_key

    def _read_params(self):
        """
        Read parameter data from the SDDS file.

        Returns
        -------
        parameters: dictionary
            Dictionary object with parameters names and values.
        """

        # if self.parameters:
        #     return self.parameters
        # else:
        #     pass

        try:
            self.parsef
        except AttributeError:
            self._read_header()
            self._parse_header()
            if self.verbose:
                print("Header data read and parsed.")

        self.parameters = {}

        # Reset pointer back to beginning of binary data to start readin there
        # self.openf.seek(self.pointer)

        param_data = ()
        for key in self.param_key:
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

    def create_paramater(self, parName, parData, parType,
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
                    # split boolean to avoid the hassle of testing existance of arrays separately
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