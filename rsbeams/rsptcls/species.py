from scipy.constants import c


class Species:
    """
    For holding species data in rsbeams.rsdata.switchyard.Switchyard
    """
    def __init__(self, coordinates, charge=None, mass=None, total_charge=None):
        self.x = coordinates[:, 0]
        self.ux = coordinates[:, 1]
        self.y = coordinates[:, 2]
        self.uy = coordinates[:, 3]
        self.ct = coordinates[:, 4]
        self.pt = coordinates[:, 5]
        self.charge = charge
        self.mass = mass
        self.total_charge = total_charge
