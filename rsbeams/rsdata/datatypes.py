from .utils import data_types, _shlex_split
import numpy as np


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


# Available namelists from SDDS standard
supported_namelists = {'&parameter': Parameter, '&column': Column, '&data': Data, '&array': Array, '&associate': Associate}