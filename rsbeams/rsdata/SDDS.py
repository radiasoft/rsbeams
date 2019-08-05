from future.builtins import str
import numpy as np
from struct import pack, unpack, calcsize, error
from sys import byteorder
from copy import copy


class readSDDS:
    """
    Class for reading SDDS data files.
    Usage:
        After initialization parameter data can be read out using the 'read_params' method
        and column data through the 'read_columns' method.

    Caveats:
        - System is assumed little-endian
        - Data stored in binary
        - No array data (only parameters and columns)
        - Strings stored as fixed are not included in parameters (strings stored in binary data are included)
        - Files that store string data in columns are not currently supported
        - Multipage SDDS files are nt currently supported
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
        # TODO: the storage of columns and parameters as attributes is inconsistent. THey should probably just both be local variables, I don't think they need to be accessed excpet for debugging purposes.
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
                    parameter_position +=2
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
                self.param_names.append(param[i0:i0+ie])

        self.param_size = calcsize(self.param_key[0])
        self.column_size = calcsize(self.column_key)

        if self.verbose:
            print("Parameter unpack size: %s bytes \nColumn unpack size: %s bytes".format((self.param_size, self.column_size)))

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

        for i in range(self.row_count):
            self.columns.append(unpack(self.column_key, self.openf.read(self.column_size)))

        return np.asarray(self.columns)

    def _get_row_count(self, here=False):
        """
        Get row count on a page.
        Args:
            here (boolean): If False then get the row count of the first page. Else try to read at current position.

        Returns:

        """

        if not here:
            self.openf.seek(self.header_end_pointer)
        self._row_count = unpack('i', self.openf.read(calcsize('i')))[0]
        if not here:
            self.openf.seek(self.header_end_pointer)

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
        self._get_row_count()
        header_size = self.header_end_pointer
        param_size = self.param_size
        # Need the size of all rows and columns (not just one row) for this
        col_size = self.column_size * self._row_count

        # Read all requested pages
        for page in pages:

            # check if pages were skipped and update position
            expected_position = page * (param_size + col_size) + header_size
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
columnSDDS = """&column name=col%d, type=double, &end\n"""
parameterSDDS = """&parameter name=col%d, type=double, &end\n"""
columnAttributeStr = ['name=', 'type=', 'units=', 'symbol=', 'format_string=', 'description=']
parameterAttributeStr = ['name=', 'type=', 'units=', 'symbol=', 'format_string=', 'description=']


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

    key_indentity = {'double': 'd', 'short': 's', 'long': 'i'}

    def __init__(self, page=1, readInFormat='numpyArray'):
        """
        Initialize SDDS object for storing parameter/column data and writing.
        """
        self.format = readInFormat
        self.endianness = byteorder
        self.page = page
        self.columns = []
        self.columnData = []
        self.columnAttributes = []
        self.parameters = []
        self.parameterData = []
        self.parameterAttributes = []
        self.column_key = '='

    def create_column(self, colName, colData, colType, colUnits='', colSymbol='', colFormatStr='', colDescription=''):
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
        """

        self.columns.append(colName)
        self.columnAttributes.append([])
        self.columnAttributes[-1].append(colName)
        self.columnAttributes[-1].append(colType)
        self.columnData.append(colData)

        try:
            self.column_key += self.key_indentity[colType]
        except KeyError:
            print("Not a Valid Data Type")

        for attribute in (colUnits, colSymbol, colFormatStr, colDescription):
            if attribute:
                self.columnAttributes[-1].append(attribute)
            else:
                self.columnAttributes[-1].append('')

    def create_param(self, parName, parData, parType, parUnits='', parSymbol='', parFormatStr='', parDescription=''):
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
        """

        self.parameters.append(parName)
        self.parameterAttributes.append([])
        self.parameterAttributes[-1].append(parName)
        self.parameterAttributes[-1].append(parType)
        self.parameterData.append(parData)
        for attribute in (parUnits, parSymbol, parFormatStr, parDescription):
            if attribute:
                self.parameterAttributes[-1].append(attribute)
            else:
                self.parameterAttributes[-1].append('')

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
        columnString = ''
        parameterString = ''
        outputFile = open(fileName, 'w')

        if len(self.columnData) > 1:
            try:
                columnDataPrint = np.column_stack(self.columnData)
            except ValueError:
                print('ERROR: All columns on a page must have same length')
                return
        elif len(self.columnData) == 1:
            columnDataPrint = self.columnData[0]
        else:
            columnDataPrint = np.empty([0])

        # Begin header writeout
        outputFile.write(self.sddsIdentifier)
        outputFile.write('!# %s-endian\n' % self.endianness)

        for parameter in self.parameterAttributes:
            for attribute in zip(parameter, parameterAttributeStr):
                if attribute[0]:
                    parameterString = ''.join((parameterString, '%s%s, ' % (attribute[1], attribute[0])))
            outputFile.write('&parameter %s&end\n' % parameterString)
            parameterString = ''

        for column in self.columnAttributes:
            for attribute in zip(column, columnAttributeStr):
                if attribute[0]:
                    columnString = ''.join((columnString,'%s%s, ' % (attribute[1], attribute[0])))
            outputFile.write('&column %s &end\n' % columnString)
            columnString = ''

        outputFile.write('&data mode=%s, &end\n' % self.dataMode)

        # Begin data writeout
        for parameter in self.parameterData:
            if self.dataMode == 'ascii':
                outputFile.write('%s\n' % parameter)
            if self.dataMode == 'binary':
                outputFile.write('%s' % (pack('d', parameter)))
                
        if self.dataMode == 'ascii':
                outputFile.write('%s\n' % columnDataPrint.shape[0])
        if self.dataMode == 'binary':
                outputFile.write('%s' % pack('I', columnDataPrint.shape[0]))

        if self.dataMode == 'ascii':
            np.savetxt(outputFile, columnDataPrint)
        elif self.dataMode == 'binary':
            for row in columnDataPrint:
                outputFile.write('%s' % (pack(self.column_key, *row)))
        # elif self.dataMode == 'binary':
        #     columnDataPrint.astype('float64').tofile(outputFile)
        else:
            print("NOT A DEFINED DATA TYPE")

        outputFile.close()


if __name__ == '__main__':
    ff = readSDDS('run_setup.output.sdds')
    p, c = ff.read(pages=[0, 1, 2, 4])
    print(p, c)
