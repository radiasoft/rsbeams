import argparse as arg
import numpy as np
from scipy.constants import c, m_e, physical_constants

m_e_ev = physical_constants['electron mass energy equivalent in MeV'][0] * 1e6
m_e_kg = m_e
kg_per_ev = m_e_kg / m_e_ev

m = 2.

parser = arg.ArgumentParser()

input = parser.add_mutually_exclusive_group()
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
    def __init__(self, start_parser):
        print("started")
        args = start_parser.parse_args()
        self.args = {key: getattr(args, key) for key in vars(args)}

        self.startup = {"momentum": self.start_momentum,
                        "velocity": self.start_velocity,
                        "energy": self.start_energy,
                        "kinetic": self.start_kenergy,
                        "betagamma": self.start_betagamma,
                        "beta": self.start_beta,
                        "gamma": self.start_gamma}

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
                        "mass_unit": self.args["mass_unit"],
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

    @staticmethod
    def start_velocity(velocity):
        beta = velocity / c

        return beta, 1. / np.sqrt(1 - beta**2)

    @staticmethod
    def start_gamma(beta):
        beta = beta
        return beta, 1. / np.sqrt(1 - beta**2)

    @staticmethod
    def start_beta(gamma):
        beta = np.sqrt(1 - 1 / gamma**2)
        return beta, gamma

    @staticmethod
    def start_betagamma(betagamma):
        beta = betagamma / np.sqrt(1 + betagamma**2)

        return beta, 1. / np.sqrt(1 - beta**2)

    def start_momentum(self, momentum):
        normalized_momentum = momentum / (self.mass * c)
        beta = normalized_momentum / np.sqrt(1 + normalized_momentum**2)
        gamma = 1 / np.sqrt(1 - beta**2)

        return beta, gamma

    def start_energy(self, energy):
        gamma = energy / (self.mass * c**2)
        beta = np.sqrt(1. - 1 / gamma**2)

        return beta, gamma

    def start_kenergy(self, kenergy):
        gamma = kenergy / (self.mass * c**2) + 1.

        return np.sqrt(1. - 1 / gamma**2), gamma

    def calculate_momentum(self, beta, gamma, **kwargs):
        return beta * gamma * self.mass * c

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
    run_converter = Converter(start_parser=parser)
    run_converter()