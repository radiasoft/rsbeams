from rsbeams.rsstats import kinematic
import pytest
"""
top level cases: outputs / inputs
    for each: particle type 
     Read properties for particle type(mass, bet)
     INSTANTIATE Converter
        For key in result actually test values

"""



@pytest.fixture
def build_converter(request):
    vals = request.param
    c = kinematic.Converter(mass=mass, beta=beta, output_unit=output_unit,
                            input_unit=input_unit, mass_unit=input_mass_unit)
    return c


# @pytest.mark.parametrize("mass", [0.511e6, 938e6])
# @pytest.mark.parametrize("beta", [0.5, 0.75])
# def test_meta(mass, beta):
#     print(mass, beta)

@pytest.mark.parametrize("build_converter", [0.511e6, 938e6], indirect=True)
@pytest.mark.parametrize("beta,key,expected", [[0.5, 'beta', 0.500], [0.75,'beta', 0.7500]])
def test_converter_build(build_converter, beta, key, expected,
                         input_unit='eV', output_unit='eV', input_mass_unit='eV'):
    c = build_converter
    assert c[key] == expected
