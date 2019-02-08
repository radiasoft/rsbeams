from .Element import Element
from sympy import symbols, cosh, sinh, sqrt, lambdify
from sympy.matrices import Matrix
from ruamel import yaml


class StructuredBeamline(object):
    """
    Holds Element objects and may contain other StructuredBeamline objects.
    """
    # TODO: Change main control sequence and beamline object to allow for multiple layers of sublines.
    K1, L = symbols('K1 L')
    matrices = {
        # Transverse matrices for elements
        'quadrupole': Matrix(([cosh(sqrt(-K1) * L), sinh(sqrt(-K1) * L) / sqrt(-K1), 0, 0],
                              [sqrt(-K1) * sinh(sqrt(-K1) * L), cosh(sqrt(-K1) * L), 0, 0],
                              [0, 0, cosh(sqrt(K1) * L), sinh(sqrt(K1) * L) / sqrt(K1)],
                              [0, 0, sqrt(K1) * sinh(sqrt(K1) * L), cosh(sqrt(K1) * L)])),
        'drift': Matrix(([1, L, 0, 0],
                         [0, 1, 0, 0],
                         [0, 0, 1, L],
                         [0, 0, 0, 1]))
    }
    elements = {}

    def __init__(self):
        self.beamline_name = None
        self.sequence = []
        self._length = self._get_length

    @property
    def length(self):
        return self._length()
    @length.setter
    def length(self, *arg, **kwargs):
        raise AttributeError("You cannot change beamline length")

    def set_beamline_name(self, beamline_name):
        self.beamline_name = beamline_name

    def add_element(self, element_name, element_type, element_parameters, subline=False):
        """
        Insert a new element at the end of a beamline.
        Can automatically be added to a subline if subline=True.
        If `element_name` already exists a line is created to it instead. `element_type` and `element_parameters`
        will be ignored in this case.
        Args:
            element_name: (str) Name of the element.
            element_type: (str) Type of element. Ignored if element is a copy.
            element_parameters: (dict) Parameters to give to the element. Ignored if element is a copy
            subline: (bool) Should the element be appended to a subline at the end of the current beamline?

        Returns:

        """
        if element_name in self.elements:
            print('Element {} already exists. Inserting copy.'.format(element_name))
            the_element = self.elements[element_name]
        else:
            the_element = Element(element_name, element_type, **element_parameters)
            self.elements[element_name] = the_element
        if subline:
            assert isinstance(self.sequence[-1], StructuredBeamline), "Last element is not type StructuredBeamline \
                                                                      subline must be false"
            self.sequence[-1].sequence.append(the_element)
        else:
            self.sequence.append(the_element)

    def add_beamline(self, name):
        self.sequence.append(StructuredBeamline())
        self.sequence[-1].set_beamline_name(beamline_name=name)

    def save_beamline(self, filename):
        """
        WARNING: This will only work correctly on Python 3

        Create a YAML file of the beamline.
        Warning: This process does not preserve sub-line structure
        :param filename:
        :return:
        """
        # Return elements preserving sub-line structures
        def return_elements(line):
            the_beamline = []
            for ele in line.sequence:
                if isinstance(ele, Element):
                    the_beamline.append({'name': ele.name, **ele.parameters, 'type': ele.type})
                else:
                    the_beamline.append(return_elements(ele))
            return the_beamline

        beamline = return_elements(self)
        with open(filename, 'w') as outputfile:
            yaml.dump(beamline, outputfile, default_flow_style=False, Dumper=yaml.RoundTripDumper)

    def load_beamline(self, filename):
        if len(self.sequence) != 0:
            print("Cannot load a new beamline.\nThis StructuredBeamline is not empty.")
            return
        elements = yaml.load(open(filename, 'r'), Loader=yaml.Loader)

        def create_beamline(elements, beamline):
            for element in elements:
                if isinstance(element, dict):
                    beamline.add_element(element['name'], element['type'],
                                     {k: v for k, v in element.items()
                                      if (k != 'type') and (k != 'name')})
                else:
                    self.add_beamline(name=None)
                    create_beamline(element, self.sequence[-1])

        create_beamline(elements, self)

    def get_beamline_elements(self):
        """
        Returns a generator object containing all elements, in order, from the beamline and any sub-beamlines
        it contains.
        :return:
        """

        def generate_beamline(element):
            if isinstance(element, (StructuredBeamline, list)):
                try:
                    element = element.sequence
                except AttributeError:
                    pass
                for value in element:
                    for subvalue in generate_beamline(value):
                        yield subvalue
            else:
                yield element

        return generate_beamline(self.sequence)

    def _get_length(self):
        length = 0.0
        for ele in self.get_beamline_elements():
            try:
                length += ele.parameters['L']
            except KeyError:
                pass

        return length

    def edit_element(self, element, parameter, value, warnings=True):
        # TODO: Assumes all elements have unique names or that you want to edit all elements of the same name.
        """
        Change one or multiple parameters of an element.

        :param element: (int or str) Position in the beamline or name of the element to change.
        :param parameter: (str or list) Name of names (as list of strings) of the parameters to change.
        :param value: If a list was given for parameter must be of equal length. Otherwise it is left to the user
        to ensure that the appropriate type of value is assigned here.
        :param warnings: (bool) Print alert if a new parameter is created.
        :return: None
        """
        assert type(element) == str or type(element) == int, "element must be a string or int"

        if type(parameter) != list:
            parameter = [parameter, ]
        if type(value) != list:
            value = [value, ]

        if type(element) == str:
            eles = [i for i, ele in enumerate(self.sequence) if ele.name == element]
        elif type(element) == int:
            eles = [element]

        for index in eles:
            for p, v in zip(parameter, value):
                try:
                    self.sequence[index].parameters[p]
                    self.sequence[index].parameters[p] = v
                except KeyError:
                    self.sequence[index].parameters[p] = v
                    if warnings:
                        print("WARNING: Creating a new parameter {} for element {}".format(p, self.sequence[index].name))

    def generate_matrix(self, concatenate=True):
        elements = []
        variables = {}
        for ele in self.sequence:
            elements.append(self.matrices[ele.type])

            # If parameter has numeric value use it, otherwise prepare to lambdify
            for key, val in ele.parameters.items():
                if type(val) == str:
                    variables[val] = symbols(val)
                    elements[-1] = elements[-1].subs(key, val)
                else:
                    elements[-1] = elements[-1].subs(key, val)

        if concatenate:
            matrix = elements[-1]
            for ele in elements[-2::-1]:
                matrix *= ele
            elements = matrix

        if len(variables) > 0:
            eval_matrix = lambdify([val for val in variables.values()], elements)
            return eval_matrix, elements
        else:
            return elements