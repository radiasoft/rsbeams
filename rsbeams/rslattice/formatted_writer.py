import numpy as np
from scipy.constants import c, e, physical_constants, m_e
m_e_mev = physical_constants['electron mass energy equivalent in MeV'][0]
from .Element import defined_elements
# TODO: May ultimately need some super-type that classifies elements that normally have common handling
#     For instance cavities and bending magnets have particular special handling, Quadrupoles to a much lesser extent.


class Writer:

    def __init__(self, beamline):
        self.beamline = beamline

        self.cavity_preference = 'rfca'

        self.E0 = 0.
        self.centralE = 0.  # Total Energy
        self.particle_mass = m_e_mev
        self.input_magnet_strength = 'normalized'
        self.cavities_change_energy = True
        self.cavity_waveform = np.sin
        self.N = 2

    def cavity_handler(self, element):
        if self.cavities_change_energy:
            # TODO: There should probably be a factor of N in here to get correct energy increase
            deltaE = self.N * element.parameters['volt'] * self.cavity_waveform(element.parameters['phase'] * np.pi / 180.)

            self.centralE += deltaE * 1e-6  # convert deltaE to MV
            print("New central Energy after cavity: {} MeV".format(self.centralE))

    def get_brho(self):
        e_energy = self.centralE  # MeV
        gamma = e_energy / self.particle_mass
        beta = np.sqrt(1 - 1 / gamma ** 2)
        e_momentum = beta * gamma * (1e6 * self.particle_mass * e / c**2) * c  # beta * gamma * mass[kg] * c[m/s]
        b_rho = e_momentum / (self.N * e)  # kg*m/s / C  == T*m

        return b_rho

    def quadrupole_handler(self, element):
        normalized_k1 = element.parameters['k1'] / self.get_brho()

        return normalized_k1


class ElegantWriter(Writer):
    element_template = "{name}: {type},"
    parameter_template = " {parameter}={value},"
    line_template = "{line_name}: LINE=("

    def __init__(self, beamline):
        super(ElegantWriter, self).__init__(beamline)

    def write_file(self, output_format, filename, E0=None, input_magnet_strength='normalized', prefix=''):
        """

        Args:
            output_format: Output code to format for. Options currently:
                elegant
            E0: Starting energy if required.
                Required if converting between `magnet_strength` representations.
            input_magnet_strength: How magnet strengths are defined in the BeamLine representation. May be:
                'normalized': gradient / Brho for quadrupole
                'absolute': gradient in T/m

        Returns:

        """
        self.prefix = prefix
        assert output_format == 'elegant', "elegant is only output format recognized currently"
        assert input_magnet_strength == 'normalized' or input_magnet_strength == 'absolute', "magnet_strength no accepted"
        if format == 'elegant' and input_magnet_strength != 'normalized':
            assert E0 is not None, "E0 must be defined"

        if E0:
            self.E0 = E0
            self.centralE += self.E0
        print("Starting central Energy: {} MeV".format(self.centralE))

        with open(filename, 'w') as ff:
            for ele in self.beamline.get_beamline_elements():
                if ele.type == 'quadrupole':
                    new_K1 = self.quadrupole_handler(ele)
                    params = ele.parameters.copy()
                    params['k1'] = new_K1
                    new_line = self.element_format(ele, alternate_parameters=params) + '\n'
                    ff.write(new_line)
                elif ele.type == 'rfca':
                    self.cavity_handler(ele)
                    new_line = self.element_format(ele) + '\n'
                    ff.write(new_line)
                else:
                    new_line = self.element_format(ele) + '\n'
                    ff.write(new_line)

            ff.write(self.line_format())

    def element_format(self, element, alternate_parameters=None):
        # TODO: SHould reference element parameter templates and warn when a templated parameter is not filled
        if alternate_parameters:
            element_parameters = alternate_parameters.copy()
        else:
            element_parameters = element.parameters.copy()

        for par in element_parameters.copy():
            if par not in defined_elements[element.type]:
                element_parameters.pop(par)

        name = self.element_template.format(name=self.prefix + element.name, type=element.type)
        parameters = [self.parameter_template.format(parameter=p, value=v) for p, v in element_parameters.items()]
        element_string = name + ''.join(parameters)

        return element_string[:-1]

    def line_format(self):
        line = ''
        for ele in self.beamline.get_beamline_elements():
            line += ',' + self.prefix +str(ele.name)

        line_str = self.line_template.format(line_name=self.beamline.name) + line[1:] + ')'

        return line_str









