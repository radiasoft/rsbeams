from .Beamline import StructuredBeamline

"""
Convert elegant lte beamlines to the floor coordinate system of OPAL.
Translate element types to OPAL from elegant syntax.
"""

# TODO: Add a beamline length parameter
# TODO: Add a copy function

import numpy as np
digits = '0123456789'

"""
Notes on sanitizing input file and element search assumptions:

lte files have 4 types of possible entries on a line:
1. Element definition
2. Beamline definition
3. Comment
4. RPN definition

Delete comment lines and comments at end of lines
Concatenate lines with continuation characters (&)
Fine beamline definitions and remove to new list
Only element definitions left in lattice_definition
Should be safe to search for exact element name (with quotation marks, if any, stripped) within lattice_definition
"""


class BeamlineParser(object):
    """
    Class that will parse a .lte file and return a StructuredBeamline object.
    May be used to then construct an equivalent beamline in other input formats.

    Parser assumes that the .lte file can be read by elegant without error.

    Missing features:
     '-' operator for reversing beamline direction
    """

    def __init__(self, filename, beamline_name):
        self.filename = filename
        self.beamline_name = beamline_name.lower()
        self.beamline = StructuredBeamline(beamline_name)  # Holds StructuredBeamline object
        self.beamline_string = ''
        self.lines = {}
        self.rpn_variables = {}

        try:
            with open(filename, 'r') as open_file:
                self.lattice_definition = open_file.readlines()
        except IOError:
            print("File could not be read")

        self._sanitize_lattice_definition()
        self._replace_stored_variables()

    def __call__(self):
        """
        Runs parser for set beamline and returns StructuredBeamline object
        :return:
        """
        self.find_beamlines()
        self.parse_set_beamline()

        return self.beamline

    def set_beamline(self, beamline_name):
        """
        Set name of beamline to be parsed.
        :param beamline_name:
        :return:
        """
        self.beamline_name = beamline_name

    def find_beamlines(self):
        """
        Find all beamline definitions in file.
        :return:
        """
        if not self.beamline_name:
            return "Set Beamline"

        for linenumber in range(len(self.lattice_definition) - 1, -1, -1):
            if self.lattice_definition[linenumber].find('line=(') > -1:
                def_start = self.lattice_definition[linenumber].find('=') + 1
                self.lines[self._read_name(self.lattice_definition[linenumber])] = \
                    self.lattice_definition[linenumber][def_start:]

                if self._read_name(self.lattice_definition[linenumber]) == self.beamline_name:
                    def_start = self.lattice_definition[linenumber].find('=') + 1
                    self.beamline_string = self.lattice_definition[linenumber][def_start:]

                self.lattice_definition.pop(linenumber)

    def parse_set_beamline(self):
        """
        Parses the beamline specified by the beamline_name attribute and creates a StructuredBeamline object
        held by the beamline attribute.
        :return:
        """
        assert self.beamline_name in self.lines, "The requested beamline is not found"

        beamline_string = self.beamline_string[1:]  # Consume opening '('
        coefficient = ''
        lines_list = []  # Holds beamline definitions if recursion is necessary
        flag_subline = False  # Set if recursion begins

        """
        Order of operations for parsing:
            - Determine if element is repeated
            - If Yes, then expand element and replace with explicit list
            - Find element name
            - Determine if element is another beamline or a single element
            - If beamline, then store current beamline and restart
            - If a single element, then find element definition in file
            - Extract element parameters
            - Remove element name from front of string
            - Check if string is empty
            - If current beamline string not empty then restart process
            - If empty check if an beamlines exist in storage
            - If beamline then set it as active and restart
            - If no beamline exists in storage and current is empty then exit
        """

        while beamline_string.strip():
            # Check for leading coefficient on element
            while beamline_string[0] in digits:
                coefficient += beamline_string[0]
                beamline_string = beamline_string[1:]
            # All digits of coefficient will be consumed

            if coefficient:
                beamline_string = self._expand_line(beamline_string, int(coefficient))
                coefficient = ''
            if beamline_string.find(',') > -1:
                next_element_end = beamline_string.find(',')
            elif beamline_string.find(')') > -1:
                next_element_end = beamline_string.find(')')

            element_name = beamline_string[:next_element_end].strip()
            # If element is a line definition then start recursion process
            if element_name in self.lines:
                print("Starting on line:", element_name)
                if beamline_string.find(',') > -1:
                    beamline_string = beamline_string[
                                      beamline_string.find(',') + 1:]  # Consume element name and ','
                elif beamline_string.find(')') > -1:
                    beamline_string = beamline_string[
                                      beamline_string.find(')') + 1:]  # Consume element name and ')'
                if beamline_string.strip():
                    lines_list.append(beamline_string)

                # Set beamline_string to be new beamline
                beamline_string = self.lines[element_name]
                beamline_string = beamline_string[1:]  # Consume opening '('
                self.beamline.add_beamline(name=element_name)
                flag_subline = True
                continue

            # Strip out quotations in case someone wasn't consistent about using them around the element name...
            element_name = element_name.replace('\'', '')
            element_name = element_name.replace('"', '')

            # Find definition of element in lattice definition file
            element_definition = None
            for line in self.lattice_definition:
                line = line.replace('"', '')
                line = line.replace('\'', '')
                if line.lstrip().find('{}:'.format(element_name)) == 0:
                    element_definition = line

            # Add element to appropriate StructuredBeamline
            try:
                element_type, element_parameters = self._find_parameters(element_definition)
            except (NameError, TypeError):
                print("\n\nElement Definition for {} not found".format(repr(element_name)))
                raise

            self.beamline.add_element(element_name, element_type, element_parameters, subline=flag_subline)
            beamline_string = beamline_string[next_element_end + 1:]  # Consume element name

            # Catch dangling bracket
            if beamline_string.rstrip() == ')':
                beamline_string = ''

            # Move up 1 level in recursion
            if not beamline_string.strip():
                if len(lines_list) != 0:
                    beamline_string = lines_list[-1]
                    lines_list.pop()
                    flag_subline = False

    def _expand_line(self, beamline_string, coefficient):
        """
        Replaces successive elements of form N*(element/line list) or N*element or line with
        explicit members.
        :param beamline_string:
        :param coefficient:
        :return:
        """
        beamline_string = beamline_string[1:]  # Consume '*'
        subline = ''

        # Case 1: Starts with (
        if beamline_string[0] == '(':
            beamline_string = beamline_string[1:]  # Consume '('
            subline_end = beamline_string.find(')')

            for _ in range(coefficient):
                subline += beamline_string[:subline_end]
                subline += ','
            beamline_string = subline + beamline_string[subline_end + 2:]
            # TODO: This may overwrite ) if the final element, but I may not care for the purpose of the parser
        # Case 2: Not enclosed by ()
        else:

            # Catch if we are at the end of the enclosing beamline definition
            subline_end = beamline_string.find(',')
            if subline_end > -1:
                pass
            else:
                subline_end = beamline_string.find(')')

            for _ in range(coefficient):
                subline += beamline_string[:subline_end]
                subline += ','
            beamline_string = subline + beamline_string[subline_end + 1:]  # No bracket to skip, just +1

        return beamline_string

    def _find_parameters(self, beamline_string):
        """
        Finds element type and keyword parameter values
        :param beamline_string:
        :return:
        """

        type_start = len(beamline_string) - (beamline_string[::-1].find(':') + 1)
        type_end = beamline_string.find(',')

        element_type = beamline_string[type_start + 1:type_end].strip()
        beamline_string = beamline_string[type_end + 1:]

        element_parameters = {}

        # Iterate through 'keyword' parameters of the element
        while beamline_string.partition('=')[1]:
            beamline_string = beamline_string.partition('=')
            key = beamline_string[0].strip()

            try:
                element_parameters[key] = parse_rpn(beamline_string[2])
            except ValueError:
                beamline_string = beamline_string[2].partition(',')
                try:
                    element_parameters[key] = parse_rpn(beamline_string[0])
                except ValueError:
                    # TODO: This allows setting of non-numeric parameters (which do occur),
                    # TODO: but also could allow improperly parsed numbers to be added as a string with no warning.
                    element_parameters[key] = beamline_string[0]

            beamline_string = beamline_string[2]

        return element_type, element_parameters

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
                if self.lattice_definition[linenumber].lstrip()[0] == '!':
                    self.lattice_definition.pop(linenumber)
            except IndexError:
                self.lattice_definition.pop(linenumber)
                continue

            # Remove anything after a comment midline
            if self.lattice_definition[linenumber].find('!') > -1:
                self.lattice_definition[linenumber] = \
                    self.lattice_definition[linenumber][:self.lattice_definition[linenumber].find('!')]

            try:
                if self.lattice_definition[linenumber].rstrip()[-1] == '&':
                    self.lattice_definition[linenumber] = self.lattice_definition[linenumber].rstrip()[:-1] + \
                        self.lattice_definition[linenumber + 1].strip()
                    self.lattice_definition.pop(linenumber + 1)
            except IndexError:
                print()
                print("Line continuation concatenation may have failed.")
                raise

            if self.lattice_definition[linenumber].lstrip()[0] == '%':
                variable = self.lattice_definition[linenumber].split()
                self.rpn_variables['{}'.format(variable[-1])] = parse_rpn(' '.join(variable[1:-2]))

    def _replace_stored_variables(self):
        """
        Replace any variables stored in rpn_variable with their calculated value.
        """

        for linenumber in range(len(self.lattice_definition)):
            for key in self.rpn_variables:
                if self.lattice_definition[linenumber].find(key):
                    self.lattice_definition[linenumber] = \
                        self.lattice_definition[linenumber].replace(key, str(self.rpn_variables[key]))

    def _read_name(self, line):
        end = line.find(':')

        if end < 0:
            return -1

        return line[0:end]


def parse_rpn(expression):
    # TODO: Common expression: atan, pi, sqr,
    # See http://www.aps.anl.gov/Accelerator_Systems_Division/Accelerator_Operations_Physics/manuals/SDDStoolkit-6/node99.html
    # For additional capabilities that may need to be added
    """
    RPN Calculator
    :param expression:
    :return:
    """

    stack = []

    for val in expression.split(' '):
        # Sanitize string
        val = val.replace('"', '')
        val = val.replace('\'', '')
        val = val.strip()
        # TODO: Could I use a dictionary to hold operators and their code operation to get ride of if/elif
        if val in ['-', '+', '*', '/', 'mult']:
            token1 = stack.pop()
            token2 = stack.pop()
            if val == '-':
                result = token2 - token1
            elif val == '+':
                result = token2 + token1
            elif val == '*' or val == 'mult':
                result = token2 * token1
            elif val == '/':
                result = token2 / token1
            else:
                raise RuntimeError
            stack.append(result)
        elif val in ['sqrt']:
            token1 = stack.pop()
            if val == 'sqrt':
                result = np.sqrt(token1)
            else:
                raise RuntimeError
            stack.append(result)
        else:
            if val:
                stack.append(float(val))

    return stack.pop()
