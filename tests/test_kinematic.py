import pytest
import numpy
from rsbeams.rsstats import kinematic

kinematic_list = [
        "beta",
        "momentum",
    ]

@pytest.mark.parametrize(
    ("particle_mass", "kinematic_quantity"),
    [
        (938.2723e6, {"beta": 0.500, "momentum": 541711642.6722883}),
        (510998, {"beta": 0.500, "momentum": 295025.381}),
    ],
)
class TestGroup:

    @pytest.fixture
    def converter(self,  particle_mass: float, kinematic_quantity: dict, k_converter: str) -> kinematic.Converter:

        return kinematic.Converter(mass=particle_mass,
                                   mass_unit='eV',
                                   input_unit='eV',
                                   output_unit='eV',
                                   **{k_converter: kinematic_quantity[k_converter]})

    @pytest.mark.parametrize("k_converter", kinematic_list)
    @pytest.mark.parametrize("k_test", kinematic_list)
    def test_one(self, k_test,  kinematic_quantity: dict, converter: kinematic.Converter) -> None:
        assert numpy.isclose(converter(silent=True)[k_test], kinematic_quantity[k_test])
