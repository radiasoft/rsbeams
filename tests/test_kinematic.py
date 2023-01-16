import pathlib
import pytest
import numpy
from rsbeams.rsstats import kinematic
from pykern import pkjson


DATA = pkjson.load_any(open(pathlib.Path("test_resources").joinpath('kinematic_data.json')))
KINEMATICS_OPTIONS = ["beta", "momentum", "gamma", "momentum", "energy",
                  "kenergy", "brho", "velocity", "betagamma"]
UNITS_OPTIONS = ("eV", "SI")

def get_kinematic_quantities(particle_type: str, unit: str) -> dict:
    return DATA[particle_type][unit]["kinematic_quantities"]

@pytest.fixture(params=UNITS_OPTIONS)
def converter(request, particle_type: str,
              input_unit: str,  output_unit: str, k_converter: str) -> kinematic.Converter:
    mass_unit = request.param
    kinematic_quantity = get_kinematic_quantities(particle_type, unit=input_unit)

    return kinematic.Converter(mass=DATA[particle_type][mass_unit]['mass'],
                               mass_unit=mass_unit,
                               input_unit=input_unit,
                               output_unit=output_unit,
                               **{k_converter: kinematic_quantity[k_converter]})

@pytest.mark.parametrize("k_test", KINEMATICS_OPTIONS)
@pytest.mark.parametrize("k_converter", KINEMATICS_OPTIONS)
@pytest.mark.parametrize("output_unit", UNITS_OPTIONS)
@pytest.mark.parametrize("input_unit", UNITS_OPTIONS)
@pytest.mark.parametrize("particle_type", ("proton", "electron"))
def test_ev_to_ev(particle_type: str,
                  input_unit:str, output_unit:str, k_test:str, converter: kinematic.Converter) -> None:

    kinematic_quantity = get_kinematic_quantities(particle_type, unit=output_unit)
    assert numpy.isclose(converter(silent=True)[k_test], kinematic_quantity[k_test])
