from .Beamline import StructuredBeamline
import numpy as np
from re import findall

# lengths: mm
# quad strength: T/m
# bend angles: deg
# Need to handle elements that can have different number of parameters depending on mode


class BeamlineParser(object):
    """
    Class that will parse a .lte file and return a StructuredBeamline object.
    May be used to then construct an equivalent beamline in other input formats.

    """

    def __init__(self, filename, beamline_name):
        self.filename = filename
        self.beamline_name = beamline_name.lower()
        self.beamline = StructuredBeamline(beamline_name)  # Holds StructuredBeamline object
        self.beamline_string = ''
        self.lines = {}
        self.rpn_variables = {}

        self.comment_character = '!'

        try:
            with open(filename, 'r') as open_file:
                self.lattice_definition = open_file.readlines()
        except IOError:
            print("File could not be read")

        self._sanitize_lattice_definition()

    def __call__(self):
        """
        Runs parser for set beamline and returns StructuredBeamline object
        :return:
        """
        self.find_beamlines()
        self.parse_set_beamline()

        return self.beamline

    def _sanitize_lattice_definition(self):
        """
        Performs a number miscellaneous cleanup functions.
        All run here to minimize number of times the lattice definition is cycled.

        Performs:
        Remove commented and empty lines from lte file
        Concatenate lines ending in continuation character to their start
        store rpn variables
        """

        for linenumber in range(len(self.lattice_definition) - 1, -1, -1):
            self.lattice_definition[linenumber] = self.lattice_definition[linenumber].lower()

            # Remove lines starting with comment
            try:
                if self.lattice_definition[linenumber].lstrip()[0] == self.comment_character:
                    self.lattice_definition.pop(linenumber)
            except IndexError:
                self.lattice_definition.pop(linenumber)
                continue

            # Remove anything after a comment midline
            if self.lattice_definition[linenumber].find(self.comment_character) > -1:
                self.lattice_definition[linenumber] = \
                    self.lattice_definition[linenumber][:self.lattice_definition[linenumber].find(self.comment_character)]

            # Need to include this if elegant
            # try:
            #     if self.lattice_definition[linenumber].rstrip()[-1] == '&':
            #         self.lattice_definition[linenumber] = self.lattice_definition[linenumber].rstrip()[:-1] + \
            #             self.lattice_definition[linenumber + 1].strip()
            #         self.lattice_definition.pop(linenumber + 1)
            # except IndexError:
            #     print()
            #     print("Line continuation concatenation may have failed.")
            #     raise


class Trace3d(BeamlineParser):
    elements = {
        'drift': ['l'],
        'quadrupole': ['k1', 'l', 'dx', 'dy', 'of'],
        'rfca': ['etl', 'phase', 'egf', 'dwf', 'h'],
        'edge': ['beta', 'r', 'g', 'k1', 'k2'],
        'sbend': ['angle', 'r', 'index', 'vf']
    }

    classifiers = {
        1: 'drift',
        3: 'quadrupole',
        8: 'sbend',
        9: 'edge',
        10: 'rfca'
    }

    conversions = {
        # Multiply by to get BeamLine Units (mostly SI)

        # length
        'l': 1e-3,
        # dipoles
        'r': 1e-3,
        'g': 1e-3,
        'dt': np.pi / 180.,
        'angle': np.pi / 180.,
        'beta': np.pi / 180.,
        # quadrupoles
        'dx': 1e-3,
        'dy': 1e-3,
        # cavities
        'etl': 1e6

    }

    def __init__(self, filename, beamline_name):
        super(Trace3d, self).__init__(filename, beamline_name)
        self.comment_character = ';'
        self.beamline_start_position = None
        self.beamline_end_position = None

    # def __call__(self, *args, **kwargs):
    #     self._sanitize_lattice_definition()
    #     for line in self.lattice_definition:
    #         if line.find('cmt') > -1:
    #             print(self.get_element_name(line), self.get_element_type(line))

    def get_element_name(self, line):
        name_def = 'cmt('
        def_start = line.find(name_def)
        def_end = line.find(')')

        return line[def_start + len(name_def):def_end].strip()

    def get_element_type(self, line):
        def_start = line.find('nt')
        def_char_start = line[def_start:].find('=')
        element_type = findall('\d+', line[def_start + def_char_start:])[0]

        return int(element_type)

    def get_element_parameters(self, line, type):
        val_start = line.find('a(')
        val_end = line[val_start:].find('=') + val_start

        values = findall(r"[-+]?\d*\.\d+|\d+", line[val_end:])
        names = self.elements[self.classifiers[type]]

        param_dict = {}
        for param_name, param_val in zip(names, values):
            param_dict[param_name] = float(param_val)

        if len(names) != len(values):
            print("WARNING: Element {} {} had a parameter mismatch".format(self.get_element_name(line),
                                                                           self.classifiers[self.get_element_type(line)]))

        return param_dict

    def _convert_units(self, element_parameters):
        for param, conversion in self.conversions.items():
            try:
                element_parameters[param] = element_parameters[param] * conversion
            except KeyError:
                pass

    def _standardize(self):
        # This is an awkward catchall to fix non-standard parameter definitions
        for ele in self.beamline.get_beamline_elements():
            # Add arclengths to dipoles
            if ele.type == 'sbend':
                ele.parameters['l'] = np.abs(ele.parameters['r'] * ele.parameters['angle'])

    def find_beamlines(self):
        lines = []
        for i, line in enumerate(self.lattice_definition):
            if line.find('cmt') > -1:
                lines.append(i)

        self.beamline_start_position = np.min(lines)
        self.beamline_end_position = np.max(lines)

    def parse_set_beamline(self):
        for line in self.lattice_definition[self.beamline_start_position:self.beamline_end_position]:
            name = self.get_element_name(line)
            type = self.get_element_type(line)
            parameters = self.get_element_parameters(line, type)
            self._convert_units(parameters)
            self.beamline.add_element(str(name), self.classifiers[type], parameters)
            self._standardize()


