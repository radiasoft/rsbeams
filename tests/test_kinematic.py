import pathlib
import pytest
import numpy
from rsbeams.rsstats import kinematic
from pykern import pkjson


DATA = pkjson.load_any(open(pathlib.Path("test_resources").joinpath('kinematic_data.json')))
KINEMATICS_LIST = ["beta", "momentum", "gamma", "momentum", "energy",
                  "kenergy", "brho", "velocity", "betagamma"]


def get_kinematic_quantities(particle_type: str, unit: str) -> dict:
    return DATA[particle_type][unit]["kinematic_quantities"]

@pytest.fixture
def converter(particle_type: str,  k_converter: str) -> kinematic.Converter:
    mass_unit = "eV"
    kinematic_quantity = get_kinematic_quantities(particle_type, unit=mass_unit)

    return kinematic.Converter(mass=DATA[particle_type][mass_unit]['mass'],
                               mass_unit=mass_unit,
                               input_unit='eV',
                               output_unit='eV',
                               **{k_converter: kinematic_quantity[k_converter]})


@pytest.mark.parametrize("particle_type", ("proton", "electron"))
@pytest.mark.parametrize("k_converter", KINEMATICS_LIST)
@pytest.mark.parametrize("k_test", KINEMATICS_LIST)
def test_ev_to_ev(particle_type: str, k_test:str, converter: kinematic.Converter) -> None:

    kinematic_quantity = get_kinematic_quantities(particle_type, unit="eV")
    assert numpy.isclose(converter(silent=True)[k_test], kinematic_quantity[k_test])
