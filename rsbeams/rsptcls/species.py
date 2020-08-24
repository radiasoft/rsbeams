from scipy.constants import c
from rsbeams.rsstats.kinematic import Converter


class Species:
    def __init__(self, coordinates, charge=None, mass=None, total_charge=None):
        self.x = coordinates[:, 0]  # m
        self.ux = coordinates[:, 1]  # beta*gamma
        self.y = coordinates[:, 2]  # m
        self.uy = coordinates[:, 3]  # beta*gamma
        self.ct = coordinates[:, 4]  # m
        self.pt = coordinates[:, 5]  # beta*gamma
        self.charge = charge  # physical particle charge - may be float or array of floats
        self.mass = mass  # physical particle mass - may be float or array of floats
        self.total_charge = total_charge  # physical charge of all particles in bunch

    def convert_from_elegant(self):
        self.ux = self.ux * self.pt
        self.uy = self.uy * self.pt
        self.ct = self.ct * c

    def convert_from_openpmd(self):
        beta_x, gamma_x = Converter.start_beta(self.ux / (self.mass * c))
        beta_y, gamma_y = Converter.start_beta(self.uy / (self.mass * c))
        beta_z, gamma_z = Converter.start_beta(self.uz / (self.mass * c))

        self.ux = beta_x * gamma_x
        self.uy = beta_y * gamma_y
        self.uz = beta_z * gamma_z