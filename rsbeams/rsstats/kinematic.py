#!/usr/bin/env python
import sys
import argparse as arg
import numpy as np
from scipy.constants import e, c, m_e, physical_constants
from future.utils import iteritems

# TODO: Priority #1: Need to go back to the principle that user puts in whatever units they desire
# TODO: cont.  All calculations internall are cgs and then conversion is done to put in the form user requests
# TODO: Priority #2: Comprehensive set of tests for combinations of input and output units

# ERROR: kinematic -v 299788543.885 --unit SI: This returns numbers that are in eV but the units claim they're SI

# TODO: Add functionality to the Converter class for stand-alone use in scripts/notebooks
# TODO: Add output conversion to SI units
# TODO: Dynamically adjust printout to use sensible unit scale? (i.e. not just eV but MeV, GeV, etc if more readable)
# TODO: Add check when E is given to make sure E>mc**2
m_e_ev = physical_constants['electron mass energy equivalent in MeV'][0] * 1e6
m_e_kg = m_e
ev_per_kg = m_e_ev / m_e_kg

parser = arg.ArgumentParser(description="Calculate relativistic, kinematic quantities based on input of one initial"
                                        "quantity (see options below). The input may be made in appropriate "
                                        "electron volt based unit quantity or SI unit (if set). The particle type"
                                        "defaults to electron though any desired mass may set in eV/c**2 or kg "
                                        "(if appropriate flag set). ")

input = parser.add_mutually_exclusive_group(required=True)
input.add_argument("-p", "--momentum", dest="momentum", type=float, help="Input momentum value. Default unit: eV/c")
input.add_argument("-v", "--velocity", dest="velocity", type=float, help="Input velocity value. Default unit: m/s")
input.add_argument("-E", "--energy", dest="energy", type=float, help="Input velocity value. Default unit eV")
input.add_argument("-KE", "--kenergy", dest="kenergy", type=float, help="Input kinetic energy value. Default unit eV")
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
parser.add_argument("--output_unit", dest="output_unit", choices=["SI", "eV"], default="eV",
                    help="Set output unit for mass and kinetmatics. Options are:\n"
                    "'SI' for standard SI units on all outputs.\n"
                    "'eV' for the respective electron volt based unit "
                    "for all outputs.\nDefaults to 'eV'.\n")


class Converter:
    """
    Converter works by taking the input kinematic quantity and then always calculating beta and gamma;
    all other kinematic quantity calculations are then performed in terms of beta and gamma.
    """
    def __init__(self, momentum=None, velocity=None, energy=None, kenergy=None,
                 betagamma=None, beta=None, gamma=None, mass=None,
                 mass_unit='eV', input_unit='eV', output_unit='eV',
                 start_parser=None):
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
        if start_parser:
            args = start_parser.parse_args()
            self.args = {key: getattr(args, key) for key in vars(args)}
        else:
            self.args = {k: v for k, v in iteritems(locals())}
            for key in ['mass_unit', 'input_unit', 'output_unit']:
                assert self.args[key] == 'eV' or self.args[key] == 'SI', "Units must be given as SI or eV"
            declaration = 0
            for key in ['momentum', 'velocity', 'energy', 'kenergy', 'betagamma', 'beta', 'gamma']:
                if self.args[key]:
                    declaration += 1
            assert declaration == 1, "One and only one initial kinematic quantity must be provided"

        # Convert to eV for internal use
        if self.args['input_unit'] == 'SI':
            self._unit_convert(input_unit='SI', input_dict=self.args)

        # Method to call based on the kinematic quantity the user inputs
        self.startup = {"momentum": self.start_momentum,
                        "velocity": self.start_velocity,
                        "energy": self.start_energy,
                        "kenergy": self.start_kenergy,
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
                        "p_unit": "eV/c" * (self.args['output_unit'] == 'eV') + "kg * m/s" * (self.args['output_unit'] == 'SI'),
                        "e_unit": "eV" * (self.args['output_unit'] == 'eV') + "J" * (self.args['output_unit'] == 'SI'),
                        "mass": None,
                        "mass_unit": "eV/c^2" * (self.args['output_unit'] == 'eV') + "kg" * (self.args['output_unit'] == 'SI'),
                        "input": None,
                        "input_unit": None,
                        "input_type": None}

        # Match mass unit to input units. Print mass in units the user used for input though.
        if self.args['mass_unit'] == "eV":
            if self.args['mass']:
                self.mass = self.args['mass'] * (1 * (self.args['input_unit'] == 'eV') + 1 / ev_per_kg * (self.args['input_unit'] == 'SI'))
                self.outputs["mass"] = self.args['mass']
            else:
                self.mass = m_e_ev
                self.outputs["mass"] = m_e_ev
        elif self.args['mass_unit'] == "SI":
            if self.outputs['mass']:
                self.mass = self.outputs['mass'] * (1 * (self.outputs['input_unit'] == 'SI') + 1 / ev_per_kg * (self.outputs['input_unit'] == 'eV'))
                self.outputs["mass"] = self.outputs['mass']
            else:
                self.mass = m_e_ev
                self.outputs["mass"] = m_e_ev / ev_per_kg

    def __call__(self, *args, **kwargs):
        # Set beta and gamma based on the input value
        for key, value in self.args.items():
            if value and key in self.startup:
                self.outputs["input"] = value
                self.outputs["input_type"] = key
                self.outputs["beta"], self.outputs["gamma"] = self.startup[key](value)

        # Find correct unit names for the printout
        if self.args["input_unit"] == "eV":
            input_unit = "eV" * (self.outputs["input_type"] == 'energy') + \
                         "eV" * (self.outputs["input_type"] == 'kenergy') + \
                         "eV/c" * (self.outputs["input_type"] == 'momentum') + \
                         "m/s" * (self.outputs["input_type"] == 'velocity') + ""
        else:
            input_unit = "J" * (self.outputs["input_type"] == 'energy') + \
                         "J" * (self.outputs["input_type"] == 'kenergy') + \
                         "kg*m/s" * (self.outputs["input_type"] == 'momentum') + \
                         "" * (self.outputs["input_type"] == 'beta') + \
                         "m/s" * (self.outputs["input_type"] == 'velocity') + ""
        self.outputs["input_unit"] = input_unit

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
        if self.args['output_unit'] == 'SI':
            self._unit_convert(input_unit='eV', input_dict=self.outputs)

        print(print_string.format(**self.outputs))

        return self.outputs

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
        normalized_momentum = momentum / self.mass
        beta = normalized_momentum / np.sqrt(1 + normalized_momentum**2)
        gamma = 1 / np.sqrt(1 - beta**2)

        return beta, gamma

    def start_energy(self, energy):

        gamma = energy / self.mass
        beta = np.sqrt(1. - 1 / gamma**2)

        return beta, gamma

    def start_kenergy(self, kenergy):

        gamma = kenergy / self.mass + 1.

        return np.sqrt(1. - 1 / gamma**2), gamma

    # All calculate methods are called to get necessary kinematic quantities
    def calculate_momentum(self, beta, gamma, **kwargs):

        return beta * gamma * self.mass

    def calculate_energy(self, gamma, **kwargs):

        return gamma * self.mass

    def calculate_kenergy(self, gamma, **kwargs):

        return self.mass * (gamma - 1)

    def calcuate_betagamma(self, beta, gamma, **kwargs):
        return beta * gamma

    def calculate_velocity(self, beta, **kwargs):
        return beta * c

    def _unit_convert(self, input_unit, input_dict):
        # momentum energy kenergy
        conversion_factors = {
            'momentum': 1 / e * c,
            'energy': 1 / e,
            'kenergy': 1 / e,
            'mass': 1 / e * c**2
        }
        for key in conversion_factors:
            print(key, input_dict[key], (conversion_factors[key]), 2 * (input == 'SI') - 1)
            if input_dict[key]:
                input_dict[key] = input_dict[key] * (conversion_factors[key])**(2 * (input_unit == 'SI') - 1)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser.print_help()
    run_converter = Converter(start_parser=parser)
    run_converter()
