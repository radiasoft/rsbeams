from future.builtins import str
import numpy as np
import shlex
from struct import pack, unpack, calcsize, error
from sys import byteorder
from copy import copy
from types import GeneratorType

# TODO: Would be nice to refactor the old camel case convention variables
# TODO: Add multipage write support - mostly means defining how they are input
# TODO: There may be initial nuance with the row count parameter. See no_row_counts in &data command from standard.
# TODO: Need to support additional_header_lines option (never actually seen this used though)
# TODO: Will need to test buffer read after full implementation. looks like my first tests may have led me astray.
#      multipage ascii parameter read is much slower. probably due to the need to catch comments during buffer creation
# TODO: Need to add support for SDDS versions 2-4
# TODO: Handle when readSDDS.read is called multiple times for same instance. Current behavior appends to .parameters and .columns
#       probably want to have it do nothing and print a statement to set a flag to overwrite.


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
# Types with data outside the header should be appended to the end of the list for _initialize_data_arrays
sdds_namelists = ['&associate', '&description', '&include', '&column', '&parameter', '&array']
# TODO: Find sdds use case with 'short' type
# SDDS defaults to 32 bit unsigned longs on all test systems while numpy uses 64 bits for the np.int_ in test cases
#  therefore int datatypes are hardcoded in size
data_types = {'double': np.float64, 'short': np.int16, 'long': np.int32, 'string': 'S{}', 'char': np.char}


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


class Array(Datum):
    def __init__(self, namelist):
        super().__init__(namelist)


class Data(Datum):
    def __init__(self, namelist):
        self.fields = {'mode': 'binary', 'lines_per_row': 1, 'no_row_counts': 0, 'additional_header_lines': 0}
        super().__init__(namelist)
        

class Associate(Datum):
    # This field does not appear in the official SDDS spec. It is in the SDDS source.
    #  There is reference to it in very old SDDS manuals in the ToC but no matching text
    def __init__(self, namelist):
        self.fields = {'name': '', 'filename': '', 'path': '', 'description': '', 'contents': '', 'sdds': ''}
        super().__init__(namelist)


class StructData:
    def __init__(self, data_type, max_string_length):
        self.data = None
        self.max_string_length = max_string_length
        self.data_type = data_type
        self._oldarray = False

    @property
    def data_type(self):
        return self._data_type

    @data_type.setter
    def data_type(self, data_type):
        if type(data_type[0]) == list:
            self._data_type = []
            merge_hold = [b for a in data_type for b in a]
            for descr in merge_hold:
                if descr[0] == 'record_length':
                    continue
                elif type(descr[1]) == str:
                    self._data_type.append((descr[0], 'U{}'.format(self.max_string_length)))
                else:
                    self._data_type.append(descr)
        else:
            self._data_type = data_type

    def add(self, data, extend=False):
        data_hold = None
        # Stack page data data together
        for datum in data:
            datum = self._merge(datum)
            #TODO checkif needed: self._check_type(data)
            if data_hold is None:
                data_hold = datum
            else:
                data_hold = np.concatenate([data_hold, datum])
        # If this isn't the first page then extend array down axis 1
        if self.data is None:
            # Data is reshaped for parameter setup
            # if data is fed into columns the extra dimension will be quashed
            self.data = data_hold.reshape(-1, 1)
        elif extend:
            # If 
            if len(self.data.shape) == 1:
                self.data = np.array([self.data.flatten(), data_hold])
            # Add a new page for columns
            elif self.data.shape[1] == 1:
                self.data = np.array([self.data.flatten(), data_hold])
            else:
                data_hold = data_hold.reshape(-1, 1)
                self.data = np.array(list(self.data) + [data_hold.flatten()])
        else:
            # Add a new page for rows
            self.data = np.concatenate([self.data.reshape(-1, 1),
                                        data_hold.reshape(-1, 1)], axis=0)

    def _merge(self, data):
        if len(data) == 1:
            return data[0]
        else:
            new_array = np.empty(1, dtype=self.data_type)
            for arr in data:
                names = [n for n in arr.dtype.names if n != 'record_length']
                new_array[names] = arr[names]
            return new_array

    def _check_type(self, data): #Next: Need to look at adjusting string size each time up to maximum? Probably just don't merge'
        if data.dtype != np.dtype(self.data_type):
            raise TypeError("Array Datatypes do not match")




# Available namelists from SDDS standard
supported_namelists = {'&parameter': Parameter, '&column': Column, '&data': Data, '&array': Array, '&associate': Associate}


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

    def __init__(self, input_file, verbose=False, buffer=True, max_string_length=25):
        """
        Initialize the read in.

        Parameters
        ----------
        input_file: str
            Name of binary SDDS file to read.
        verbose: Boolean
            Print additional data about detailing intermediate read in process.
        buffer: Boolean
            If true then the file is entered into memory and closed before data is read. This may result in faster
            read times in some cases but only if the file is not on the order of available system memory.
        max_string_length: Int
            Upper bound on strings that can be read in. Should be at least as large as the biggest string in the file.
            Only used for formatting read from ASCII files. For binary files the string size is determined dynamically
            during the file read.
        """
        # TODO: Clean up unused attributes at the end
        self._input_file = input_file
        self.verbose = verbose
        self.buffer = buffer
        self.openf = open(input_file, 'rb')
        self.position = 0

        self.max_string_length = max_string_length
        self.header = []
        self._variable_length_records = False
        self._data_mode = None

        self.param_key = ['=i']  # Include row count with parameters
        self.column_key = '='
        self.header_end_pointer = 0  #
        self._header_line_count = 1

        self._parameter_keys = []
        self.param_size = 0
        self.param_names = ['rowCount']
        self.parameters = None

        self._column_keys = []
        self.column_size = 0
        self.column_names = {}
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
        if _read_line(self.openf).find('SDDS1') < 0:
            # First line must identify as SDDS file
            raise Exception("Header cannot be read")
        # TODO: Need a catch for &include here to at least one level of nesting
        self._header_line_count = 1
        while True:
            namelist = []
            new_line = _read_line(self.openf)
            self._header_line_count += 1
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
                print(namelist_type, namelist, 'not parsed')  # TODO: TEMP printout
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
        # TODO: the row count will must always be created in binary mode, even if there are no other parameters
        for name in sdds_namelists[3:]:  # TODO: I have using the order of the list to do this. Just test if the proper method is there.
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
            self._parameter_keys.insert(0, [('row_counts', data_types['long'])])
        elif not self.data['&data'][0].fields['no_row_counts'] and len(self.data['&column']) > 0:
            # ASCII: count may not be included and will be at the end of the parameters
            self._parameter_keys.append([('row_counts', data_types['long'])])
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

    def read2(self, pages=None):
        # Always start after the header
        if not self.buffer:
            self.openf.seek(self.header_end_pointer)

        # Select pages to be stored during the read - if None then store everything
        # pages: internal counter for page numbers
        # user_pages: store what pages are returned 
        if pages:
            user_pages = pages
        else:
            user_pages = iter_always()
        pages = iter_always()

        # Buffer reading has no pointer, you always start at the beginning so we need to move the start offset
        if self._data_mode == 'binary':
            position = 0 + self.header_end_pointer * self.buffer
        else:
            position = 0 + self._header_line_count * self.buffer

        for page in pages:
            if not isinstance(user_pages, GeneratorType) and page > np.max(user_pages):
                print(page, user_pages)
                print('stopping here')
                break

            if self._check_file_end(position):
                # TODO: Get the logic right here. Should not trigger unless user requests a bad page.
                if page in user_pages:
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
                row_count = self.parameters['row_counts'][-1][-1]
            else:
                row_count = self._get_ascii_row_count(position)

            if row_count is 0:
                continue
            if (isinstance(user_pages, GeneratorType) or (page in user_pages)) or self._variable_length_records:
                column_data, position = self._get_column_data(self._column_keys, position, row_count)
                self._columns.add(column_data, extend=True)
            else:
                # still need to update position what would have been read
                position += np.dtype(self._column_keys[0]).itemsize * row_count

    def _get_parameter_data(self, data_keys, position):
        data_arrays = [[]]
        # TODO: The variable record lengths portion has not actually been tested yet
        for dk in data_keys:
            if self._variable_length_records:
                try:
                    record_length = data_arrays[1][-1]['record_length']
                    dk = (dk[0][0], dk[0][1].format(record_length[0]))
                except (ValueError, IndexError):
                    pass
            if self.data['&data'][0].fields['mode'] == 'ascii':
                new_array = self._get_reader()(self.openf, skip_header=position, dtype=dk, max_rows=1,
                                               comments='!', deletechars='', unpack=True)
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
                    position += 1 * (row_count + 1 * int(self.data['&data'][0].fields['no_row_counts']))
            else:
                new_array = self._get_reader()(self.openf, dtype=dk, count=row_count, offset=position)
                if self.buffer:
                    position += np.dtype(dk).itemsize * row_count
            data_arrays[-1].append(new_array)

        return data_arrays, position

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