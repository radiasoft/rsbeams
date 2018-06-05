import sys
import argparse as arg
import numpy as np
from scipy.constants import c, m_e, physical_constants

# TODO: Add functionality to the Converter class for stand-alone use in scripts/notebooks

m_e_ev = physical_constants['electron mass energy equivalent in MeV'][0] * 1e6
m_e_kg = m_e
kg_per_ev = m_e_kg / m_e_ev

m = 2.

parser = arg.ArgumentParser(description="Calculate relativistic, kinematic quantities based on input of one initial"
                                        "quantity (see options below). The input may be made in appropriate "
                                        "electron volt based unit quantity or SI unit (if set). The particle type"
                                        "defaults to electron though any desired mass may set in eV/c**2 or kg "
                                        "(if appropriate flag set). ")

input = parser.add_mutually_exclusive_group(required=True)
input.add_argument("-p", "--momentum", dest="momentum", type=float, help="Input momentum value. Default unit: eV/c")
input.add_argument("-v", "--velocity", dest="velocity", type=float, help="Input velocity value. Default unit: m/s")
input.add_argument("-E", "--energy", dest="energy", type=float, help="Input velocity value. Default unit eV")
input.add_argument("-KE", "--kenergy", dest="kinetic", type=float, help="Input kinetic energy value. Default unit eV")
input.add_argument("-bg", "--betagamma", dest="betagamma", type=float, help="Input beta*gamma value. Default unit none")
input.add_argument("-b", "--beta", dest="beta", type=float, help="Input beta value. Default unit none")
input.add_argument("-g", "--gamma", dest="gamma", type=float, help="Input gamma value. Default unit none")

parser.add_argument("-m", "--mass", dest="mass", type=float, help="Optional: Value of particle mass to use.")

parser.add_argument("--mass_unit", dest="mass_unit", choices=["SI", "eV"], default="eV",
                    help="Set mass units. Options are:\n"
                    "'SI' for standard SI units on all inputs.\n"
                    "'eV' for the respective electron volt based unit "
                    "for all inputs.\nDefaults to 'eV'.\n")
parser.add_argument("--unit", dest="input_unit", choices=["SI", "eV"], default="eV",
                    help="Set input units. Options are:\n"
                    "'SI' for standard SI units on all inputs.\n"
                    "'eV' for the respective electron volt based unit "
                    "for all inputs.\nDefaults to 'eV'.\n")


class Converter:
    """
    Converter works by taking the input kinematic quantity and then always calculating beta and gamma;
    all other kinematic quantity calculations are then performed in terms of beta and gamma.
    """
    def __init__(self, start_parser):
        """
        Class that takes in a single kinematic quantity and particle mass and returns a list of other kinematic
        quantities. Options for input and the output are:
            - Velocity
            - Beta
            - Gamma
            - Momentum
            - Normalized momentum (beta * gamma)
            - Energy
            - Kinetic Energy

        Currently only works through passing in a parser object with appropriate settings
        Args:
            start_parser: Parser object containing necessary settings.

        Returns:
            None
            Prints results
        """
        print("started")
        args = start_parser.parse_args()
        self.args = {key: getattr(args, key) for key in vars(args)}

        # Method to call based on the kinematic quantity the user inputs
        self.startup = {"momentum": self.start_momentum,
                        "velocity": self.start_velocity,
                        "energy": self.start_energy,
                        "kinetic": self.start_kenergy,
                        "betagamma": self.start_betagamma,
                        "beta": self.start_beta,
                        "gamma": self.start_gamma}

        # Store quantities needed for output or method to calculate that quantity
        self.outputs = {"momentum": self.calculate_momentum,
                        "velocity": self.calculate_velocity,
                        "energy": self.calculate_energy,
                        "kenergy": self.calculate_kenergy,
                        "betagamma": self.calcuate_betagamma,
                        "beta": None,
                        "gamma": None,
                        "p_unit": "eV/c" * (args.mass_unit == 'eV') + "kg * m/s" * (args.mass_unit == 'SI'),
                        "e_unit": "eV" * (args.mass_unit == 'eV') + "J" * (args.mass_unit == 'SI'),
                        "mass": None,
                        "mass_unit": "eV/c^2" * (args.mass_unit == 'eV') + "kg" * (args.mass_unit == 'SI'),
                        "input": None,
                        "input_unit": self.args["input_unit"],
                        "input_type": None}

        # Always use eV for internal calculations. Default to electron mass if not set by user.
        if args.mass_unit == "eV":
            if args.mass:
                self.mass = args.mass
            else:
                self.mass = m_e_ev
        elif args.mass_unit == "SI":
            if args.mass:
                self.mass = args.mass * kg_per_ev
            else:
                self.mass = m_e_ev

        self.outputs["mass"] = self.mass

    def __call__(self, *args, **kwargs):
        # Set beta and gamma based on the input value
        for key, value in self.args.items():
            if value and key in self.startup:
                self.outputs["input"] = value
                self.outputs["input_type"] = key
                self.outputs["beta"], self.outputs["gamma"] = self.startup[key](value)

        # Set all derived kinematic quantities
        for key, value in self.outputs.items():
            if callable(value):
                self.outputs[key] = value(**self.outputs)

        print_string = """
        Based on an input {input_type} of {input} {input_unit}
        For a particle with mass: {mass} {mass_unit}
        
        velocity: {velocity} m/s
        beta: {beta}
        gamma: {gamma}
        momentum: {momentum} {p_unit}
        beta * gamma: {betagamma}
        energy: {energy} {e_unit}
        kinetic energy: {kenergy} {e_unit}
        """

        print(print_string.format(**self.outputs))

    # All start methods used to convert input kinematic quantity to beta and gamma
    @staticmethod
    def start_velocity(velocity):
        """
        Calculate beta and gamma based on velocity.
        Args:
            velocity: Particle velocity in m/s.

        Returns:
            (beta, gamma)
        """
        beta = velocity / c

        return beta, 1. / np.sqrt(1 - beta**2)

    @staticmethod
    def start_gamma(gamma):
        """
        Calculate beta and gamma based on gamma.
        Args:
            gamma: Relativistic gamma, unitless.

        Returns:
            (beta, gamma)

        """

        return np.sqrt(1. - 1 / gamma**2), gamma

    @staticmethod
    def start_beta(beta):
        """
        Calculate beta and gamma based on beta
        Args:
            beta: Relavistic beta, unitless

        Returns:
            (beta, gamma)
        """

        return beta, 1 / np.sqrt(1. - beta**2)

    @staticmethod
    def start_betagamma(betagamma):
        """
        Calculate beta and gamma based on beta * gamma
        Args:
            betagamma: Normalized momentum beta * gamma

        Returns:
            (beta, gamma)
        """
        beta = betagamma / np.sqrt(1 + betagamma**2)

        return beta, 1. / np.sqrt(1 - beta**2)

    def start_momentum(self, momentum):
        normalized_momentum = momentum / (self.mass * c)
        beta = normalized_momentum / np.sqrt(1 + normalized_momentum**2)
        gamma = 1 / np.sqrt(1 - beta**2)

        return beta, gamma

    def start_energy(self, energy):
        if self.args["mass_unit"] == "SI":
            mass_factor = c**2
        else:
            mass_factor = 1.0

        gamma = energy / (self.mass * mass_factor)
        beta = np.sqrt(1. - 1 / gamma**2)

        return beta, gamma

    def start_kenergy(self, kenergy):
        if self.args["mass_unit"] == "SI":
            mass_factor = c**2
        else:
            mass_factor = 1.0

        gamma = kenergy / (self.mass * mass_factor) + 1.

        return np.sqrt(1. - 1 / gamma**2), gamma

    # All calculate methods are called to get necessary kinematic quantities
    def calculate_momentum(self, beta, gamma, **kwargs):
        if self.args["mass_unit"] == "SI":
            mass_factor = c
        else:
            mass_factor = 1.0
        return beta * gamma * self.mass * mass_factor

    def calculate_energy(self, gamma, **kwargs):
        if self.args["mass_unit"] == "SI":
            mass_factor = c**2
        else:
            mass_factor = 1.0

        return gamma * self.mass * mass_factor

    def calculate_kenergy(self, gamma, **kwargs):
        if self.args["mass_unit"] == "SI":
            mass_factor = c ** 2
        else:
            mass_factor = 1.0

        return self.mass * mass_factor * (gamma - 1)

    def calcuate_betagamma(self, beta, gamma, **kwargs):
        return beta * gamma

    def calculate_velocity(self, beta, **kwargs):
        return beta * c

if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser.print_help()
    run_converter = Converter(start_parser=parser)
    run_converter()