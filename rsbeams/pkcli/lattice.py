import numpy as np
import sirepo.lib
from rsbeams.rslattice.templates.lattice_importer import parse_json
from rsbeams.rslattice.formatted_writer import WarpWriter
from rsbeams.rsstats.kinematic import Converter
from scipy.constants import c, e, physical_constants, m_e
m_e_mev = physical_constants['electron mass energy equivalent in MeV'][0]


def _get_gradient(k1, energy, mass, n=1):
    """
    K1 (float): Quadrupole normalized strength in m^-(n+1)
    energy (float): total energy in MeV
    mass (float): mass in MeV/c^2
    n (int) [default=1]: Order of multipole (n=1: quadrupole, n=2: sextupole, ...)
    """
    gamma = energy / mass
    beta = np.sqrt(1 - 1 / gamma ** 2)
    e_momentum = beta * gamma * (1e6 * mass * e / c**2) * c  # beta * gamma * mass[kg] * c[m/s]
    b_rho = e_momentum / (1 * e)  # kg*m/s / C  == T*m
    gradient = k1 * b_rho**n
    return gradient

def _parameter_conversions(beamline, energy, mass=m_e_mev):
    for ele in beamline.get_beamline_elements():
        if ele.type == 'kquad':
            ele.parameters['db'] = _get_gradient(float(ele.parameters['k1']), energy, mass)
        elif ele.type == 'csbend' or ele.type == 'rbend' or ele.type == 'csbend':
            ele.parameters['rc'] = float(ele.parameters['l']) / float(ele.parameters['angle'])
        elif ele.type == 'sext':
            ele.parameters['db'] = _get_gradient(float(ele.parameters['k2']), energy, mass, n=2)

# Public Functions

def elegant_to_warp(input_file_path, output_file_path, beamline_name):
#     input_file_path = 'check_elegant/runIOTA.ele'
    input_simulation_type = 'elegant'
#     beamline_name = 'iota'
#     output_file_path = 'warp_lattice.py'

    sim_dat = sirepo.lib.Importer(input_simulation_type).parse_file(input_file_path)
    for command in sim_dat.models.commands:
        if command['_type'] == 'run_setup':
            p_central_mev = float(command['p_central_mev']) * 1e6  # in eV for Converter
            central_energy = Converter(momentum=p_central_mev)()['energy']
            print(f"Found 'p_central_mev': {p_central_mev}\nCentral Energy is: {central_energy}\n")
            break
    else:
        raise ValueError("run_setup not found in elegant file import. Cannot set 'p_central_mev'")


    beamline = parse_json(sim_dat, beamline_name)
    _parameter_conversions(beamline=beamline, energy=p_central_mev)
    print('Warning: Conversion does not account for beam momentum changes (e.g. rf cavities, radiation loss, etc.) magnet strengths may be wrong if any of these are present')
    writer = WarpWriter(beamline, input_simulation_type)
    writer.write_file(output_file_path)
