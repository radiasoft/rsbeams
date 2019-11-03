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
        'drift': ['L'],
        'quadrupole': ['K1', 'L', 'DX', 'DY', 'OF'],
        'rfca': ['VOLT', 'PHASE', 'EGF', 'CHANGE_P0', 'H'],
        'edge': ['BETA', 'R', 'G', 'K1', 'K2'],
        'sbend': ['ANGLE', 'R', 'INDEX', 'VF']
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
        'L': 1e-3,
        # dipoles
        'R': 1e-3,
        'G': 1e-3,
        'dt': np.pi / 180.,
        'ANGLE': np.pi / 180.,
        'BETA': np.pi / 180.,
        # quadrupoles
        'DX': 1e-3,
        'DY': 1e-3,
        # cavities
        'VOLT': 1e6,

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
            try:
                self.handlers[self.classifiers[type]](name, parameters, self.beamline)
            except KeyError:
                self.beamline.add_element(str(name), self.classifiers[type], parameters)
            self._standardize()  # TODO: Change this to handler

    def handle_bends(self):
        bend_edges = []
        seq = self.beamline.sequence
        for i, ele in enumerate(seq):
            if ele.type == 'sbend':
                E1 = seq[i - 1].parameters['BETA']
                E2 = seq[i + 1].parameters['BETA']
                HGAP = seq[i - 1].parameters['G'] / 2.
                FINT = seq[i - 1].parameters['K1']
                self.beamline.edit_element(i, ['E1', 'E2', 'HGAP', 'FINT'], [E1, E2, HGAP, FINT])

                bend_edges.extend([seq[i - 1], seq[i + 1]])

        for edge in bend_edges:
            seq.remove(edge)

    def handle_cavities(self):
        cavities = []
        phases = []
        seq = self.beamline.sequence
        for i, ele in enumerate(seq):
            if ele.type == 'rfca':
                cavities.append(i)
                phases.append(90. - ele.parameters['phase'] )

                if True:
                    L = seq[i - 1].parameters['L'] + seq[i + 1].parameters['L']
                    self.beamline.edit_element(i - 1, ['L',], [0.0,])
                    self.beamline.edit_element(i + 1, ['L',], [0.0,])
                    self.beamline.edit_element(i, ['L', ], [L, ])

        for cav in cavities:
            self.beamline.edit_element(cav, ['phase']*len(phases), phases)


class TraceWin(Trace3d):
    elements = {
        'drift': ['L'],
        'quadrupole': ['K1', 'L', 'DX', 'DY', 'OF'],
        'sbend': ['ANGLE', 'R', 'INDEX', 'VF'],
        'dtl_cell': ['L', 'LQ1', 'LQ2', 'GC', 'B1', 'B2', 'VOLT', 'PHASE', 'R', 'P', 'BETA_S', 'T_S', 'KTS', 'K2TS']
    }

    classifiers = {
        'drift': 'drift',
        'quad': 'quadrupole',
        'dtl_cel': 'dtl_cell'
    }

    conversions = {
        # Multiply by to get BeamLine Units (mostly SI)

        # length
        'L': 1e-3,
        # quadrupoles

        # cavities
        'LQ1': 1e-3,
        'LQ2': 1e-3,
    }

    @staticmethod
    def handle_dtl_cell(name, cell, beamline):
        q1_L = cell['LQ1']
        q2_L = cell['LQ2']
        cavity_L = cell['L'] - q1_L - q2_L

        q1_K1 = cell['B1']
        q2_K1 = cell['B2']

        cavity_PHASE = cell['PHASE'] + 90.

        beamline.add_element(name + '_q1', 'quadrupole', {'k1': q1_K1, 'l': q1_L})
        beamline.add_element(name + 'cav', 'rfca', {'l': cavity_L, 'volt': cell['VOLT'], 'phase': cavity_PHASE})
        beamline.add_element(name + '_q2', 'quadrupole', {'k1': q2_K1, 'l': q2_L})

    handlers = {
        'dtl_cell': handle_dtl_cell.__func__
    }

    def __init__(self, filename, beamline_name):
        super(TraceWin, self).__init__(filename, beamline_name)
        self.comment_character = ';'
        self.beamline_start_position = None
        self.beamline_end_position = None

        self.name_count = 0  # TRACEWIN doesn't name elements

    def get_element_parameters(self, line, type):
        values = findall(r"[-+]?\d*\.\d+|\d+", line)
        names = self.elements[self.classifiers[type]]

        param_dict = {}
        for param_name, param_val in zip(names, values):
            param_dict[param_name] = float(param_val)

        if len(names) != len(values):
            print("WARNING: Element {} {} had a parameter mismatch".format(self.get_element_name(line),
                                                                           self.classifiers[self.get_element_type(line)]))

        return param_dict

    def get_element_type(self, line):
        element_name = line.split()[0]

        return element_name.lower()

    def get_element_name(self, line):
        name = 'e' + str(self.name_count)
        self.name_count += 1

        return name

    def _standardize(self):
        pass

    def find_beamlines(self):
        # Manually set beamline file line positions for now
        pass







