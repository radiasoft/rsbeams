# TODO: No longer going to check against a defined_elements list
#  Need an object that links like element types between codes when a StructuredBeamline
#  is dumped back to file for another code to interpret

# TODO: Keeping until above is implemented
# defined_elements = {
#     'drift': ['l'],
#     'edrift': ['l'],
#     'quadrupole': ['l', 'k1'],
#     'kquadrupole': ['l', 'k1'],
#     'solenoid': ['None'],
#     'kquad': ['l', 'k1'],
#     'edge': ['None'],
#     'sbend': ['l', 'angle', 'e1', 'e2', 'hgap', 'fint'],
#     'csbend': ['l', 'angle', 'e1', 'e2'],
#     'csrcsbend': ['l', 'angle', 'e1', 'e2'],
#     'csrdrift': ['L'],
#     'rfca': ['l', 'volt', 'frequency', 'phase', 'h'],
#     'ecol': ['None'],
#     'sextupole': ['None'],
#     'marker': ['None'],
#     'bumper': ['None'],
#     'maxamp': ['None'],
#     'watch': ['None'],
#     'scraper': ['None'],
#     'pfilter': ['None'],
#     'kicker': ['None'],
#     'vkicker': ['None'],
#     'hkicker': ['None'],
#     'monitor': ['None'],  # Apparently, even though elegant's manual calls these 'moni', 'hmon', and 'vmon'
#     'hmonitor': ['None'],  # internally they have the full word 'monitor' in their name...
#     'vmonitor': ['None'],
#     'twiss': ['None'],
#     'wake': ['None'],
#     'malign': ['None'],
#     'charge': ['None'],
#     'ibscatter': ['None']}


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

    def __init__(self, element_name, element_type, **kwargs):
        # Element properties
        self.name = element_name
        self.type = self._interpret_element(element_type)
        self.parameters = AgnosticDict()
        self._edge = self._get_edge

        # Metadata
        self._beamline = None  # Top level beamline if any

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
        # TODO: because repeated element instances are all just references this will fail when called for a repeated
        # TODO:  element. Could just have the attributes point to the same location and make the Element instance
        # TODO:  different for each. Need to test this approach. Means nothing else will be shared.
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
        # This used to verify against a list of elegant elements
        # We no longer check an element at load, only dump
        # So now this just performs sanitizing functions

        name = name.lower().strip()
        return name