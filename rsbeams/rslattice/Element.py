defined_elements = {
    'drift': ['L'],
    'quadrupole': ['L', 'K1'],
    'kquad': ['L', 'K1'],
    'sbend': ['L', 'ANGLE', 'E1', 'E2'],
    'csbend': ['L', 'ANGLE', 'E1', 'E2'],
    'csrcsbend': ['L', 'ANGLE', 'E1', 'E2'],
    'csrdrift': ['L'],
    'rfca': ['None'],
    'ecol': ['None'],
    'sextupole': ['None'],
    'marker': ['None'],
    'bumper': ['None'],
    'maxamp': ['None'],
    'watch': ['None'],
    'scraper': ['None'],
    'pfilter': ['None'],
    'kicker': ['None'],
    'vkick': ['None'],
    'hkick': ['None'],
    'monitor': ['None'],  # Apparently, even though elegant's manual calls these 'moni', 'hmon', and 'vmon'
    'hmonitor': ['None'],  # internally they have the full word 'monitor' in their name...
    'vmonitor': ['None'],
    'twiss': ['None'],
    'wake': ['None'],
    'malign': ['None'],
    'charge': ['None']}


class AgnosticDict(dict):
    def __init__(self, *args):
        dict.__init__(self, args)

    def __getitem__(self, key):
        if type(key) == str:
            return dict.__getitem__(self, key.lower())

    def __setitem__(self, key, val):
        if type(key) == str:
            dict.__setitem__(self, key.lower(), val)


class Element(object):
    """
    Class for holding parameters and attributes of beamline elements as common defined in particle tracking codes.
    Element parameters are not case sensitive.
    """
    element_types = defined_elements

    def __init__(self, element_name, element_type, **kwargs):
        # Element properties
        self.name = element_name
        self.type = self._interpret_element(element_type)
        self.parameters = AgnosticDict()
        self._edge = self._get_edge

        # Metadata
        self._beamline = None  # Top level beamline if any

        # Was going to maintain list of excepted parameters
        # Feature is not working and may just be a bad idea
        # This check for not just makes sure the element type exists at all
        try:
            type_dictionary = self.element_types[self.type]
        except KeyError:
            print(self.type)
            raise

        for parameter, value in kwargs.items():
            # Might be nice to check at some point
            # assert parameter in type_dictionary, "Parameter: {} is not supported \
            #  by element type: {}".format(parameter, self.type)

            self.parameters[parameter] = value

    @property
    def edge(self):
        """
        Position of the upstream edge of the element.
        The position is given relative to the start of the top-level beamline containing this Element.
        Returns: (float) position in meters.

        """
        return self._edge()
    @edge.setter
    def edge(self, *args, **kwargs):
        raise AttributeError("You cannot directly set the element edge")

    @property
    def beamline(self):
        return self._beamline
    @beamline.setter
    def beamline(self, *args, **kwargs):
        raise AttributeError("You cannot change the element's address")

    def _get_edge(self):
        length = 0
        for ele in self._beamline.get_beamline_elements():
            if ele != self:
                try:
                    length += ele.parameters['L']
                except KeyError:
                    pass
            else:
                break

        return length

    def _interpret_element(self, name):
        # TODO: I am leaving out a check for duplicate names on the theory that we are using a valid lte file
        """
        elegant seems to allow shortning of element names
        as long is it is uniquely defined.

        This method implements that functionality.
        """
        name = name.lower().strip()
        for element in self.element_types:
            if name == element[:len(name)]:
                return element

        print("COULD NOT FIND:", repr(name))