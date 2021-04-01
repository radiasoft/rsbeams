from .Beamline import StructuredBeamline
import numpy as np
from scipy.constants import c
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
        self.global_parameters = {}

        self.comment_character = '!'

        try:
            with open(filename, 'r') as open_file:
                self.lattice_definition = open_file.readlines()
        except IOError:
            print("File could not be read")


    def __call__(self):
        """
        Runs parser for set beamline and returns StructuredBeamline object
        :return:
        """
        self._sanitize_lattice_definition()
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
                getattr(self, self.handlers[type])(name, parameters, self.beamline)
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
        'drift': ['L', 'R', 'Ry', 'RX_SHIFT', 'RY_SHIFT'],
        'quadrupole': ['L', 'K1', 'R', 'THETA', 'G3U3', 'G4U4', 'G5U5', 'G6U6', 'GFR'],
        'sbend': ['ANGLE', 'R', 'INDEX', 'VF'],
        'dtl_cell': ['L', 'LQ1', 'LQ2', 'GC', 'B1', 'B2', 'VOLT', 'PHASE', 'R', 'P', 'BETA_S', 'T_S', 'KTS', 'K2TS'],
        'marker': [],
        'rfca': ['MODE', 'N', 'BETA', 'VOLT', 'PHASE', 'R', 'P', 'KE0TI', 'KE0TO', 'DZI', 'DZO'],
        'freq': ['FREQUENCY', ],
        'field_map': ['TYPE', 'L', 'PHASE', 'R', 'KB', 'KE', 'KI', 'KA', 'FILE']
    }

    classifiers = {
        'drift': 'drift',
        'quad': 'quadrupole',
        'dtl_cel': 'dtl_cell',
        'ncells': 'rfca',
        'set_sync_phase': 'marker',
        'freq': 'freq',
        'field_map': 'field_map'  # This is really dependent on field_map settings
    }

    conversions = {
        # Multiply by to get BeamLine Units (mostly SI)

        # length
        'L': 1e-3,
        # quadrupoles

        # cavities
        'LQ1': 1e-3,
        'LQ2': 1e-3,
        'frequency': 1e6
    }

    def handle_dtl_cell(self, name, element, beamline):
        q1_L = element['LQ1']
        q2_L = element['LQ2']
        cavity_L = element['L'] - q1_L - q2_L

        q1_K1 = element['B1']
        q2_K1 = element['B2']

        cavity_PHASE = element['PHASE'] + 90.

        beamline.add_element(name + '_q1', 'quadrupole', {'k1': q1_K1, 'l': q1_L})

        try:
            frequency = self.global_parameters['frequency']
        except KeyError:
            frequency = 1.0
        beamline.add_element(name + 'cav', 'rfca',
                             {'l': cavity_L, 'volt': element['VOLT'],
                              'phase': cavity_PHASE, 'frequency': frequency})
        beamline.add_element(name + '_q2', 'quadrupole', {'k1': q2_K1, 'l': q2_L})

    def handle_ncells(self, name, element, beamline):
        assert self.global_parameters['frequency'], "NCELL must have a frequency set to determine length"
        wavelength = c / self.global_parameters['frequency']
        length = element['BETA'] * wavelength * element['N']
        cavity_PHASE = element['PHASE'] + 90.
        cavity_VOLT = element['VOLT'] * length
        cavity_FREQ = self.global_parameters['frequency']
        beamline.add_element(name + 'cav', 'rfca',
                             {'l': length, 'volt': cavity_VOLT, 'phase': cavity_PHASE,
                              'frequency': cavity_FREQ})

    def handle_freq(self, name, element, beamline):
        self.global_parameters['frequency'] = element['frequency']

    def handle_field_map(self, name, element, beamline):
        # Very crude implementation of this element. We will assume just 1D RF field for now.
        # We will assume the field map file has max/min of 1.0 / -1.0.

        frequency = self.global_parameters['frequency']
        volt = element['VOLT'] * 1e6  # field map files store E-field in MV/m

        beamline.add_element(name + 'cav', 'rfca',
                             {'l': element['L'], 'volt': volt, 'phase': element['PHASE'],
                              'frequency': frequency})


    handlers = {
        'dtl_cel': 'handle_dtl_cell',
        'freq': 'handle_freq',
        'ncells': 'handle_ncells',
        'field_map': 'handle_field_map'
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







