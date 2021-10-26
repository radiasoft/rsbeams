import numpy


class Species:
    """
    For holding species data in rsbeams.rsdata.switchyard.Switchyard
    """

    def __init__(self, coordinates, charge=None, mass=None, total_charge=None):
        # since coordinates is mutable this isn't really very strong protection
        self._coordinates = coordinates
        self.charge = charge
        self.mass = mass
        self.total_charge = total_charge

    @property
    def coordinates(self):
        return self._coordinates

    @property
    def x(self):
        return self._coordinates[:, 0]

    @property
    def ux(self):
        return self._coordinates[:, 1]

    @property
    def y(self):
        return self._coordinates[:, 2]

    @property
    def uy(self):
        return self._coordinates[:, 3]

    @property
    def ct(self):
        return self._coordinates[:, 4]

    @property
    def pt(self):
        return self._coordinates[:, 5]

    @property
    def macroparticle_count(self):
        return self._coordinates.shape[0]
