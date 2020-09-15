import numpy as np
from pykern.pkcollections import PKDict
from sirepo.template.madx_converter import _build_field_map
from scipy.constants import c, e, physical_constants, m_e
m_e_mev = physical_constants['electron mass energy equivalent in MeV'][0]
# TODO: May ultimately need some super-type that classifies elements that normally have common handling
#     For instance cavities and bending magnets have particular special handling, Quadrupoles to a much lesser extent.

# Necessary elegant - warp mapping written in madx_converter model
_FIELD_MAP = PKDict(
    elegant=[
        ['Drft',
            ['DRIF', 'l'],
            ['CSRDRIFT', 'l'],
            ['EDRIFT', 'l'],
            ['LSCDRIFT', 'l'],
        ],
        ['Bend',
            ['CSBEND', 'l', 'rc'],
            ['SBEN', 'l', 'rc'],
            ['CSRCSBEND', 'l', 'rc'],
            ['KSBEND', 'l', 'rc'],
            ['NIBEND', 'l', 'rc'],
        ],
        ['Bend',
            ['RBEN', 'l', 'rc'],
            ['TUBEND', 'l', 'rc'],
        ],
        ['Quad',
            ['QUAD', 'l', 'db'],
            ['KQUAD', 'l', 'db'],
        ],
        ['Sext',
            ['SEXT', 'l', 'db'],
            ['KSEXT', 'l', 'db'],
        ]

        ]
)

# Significant departure from previous writer structure. Only handles output file formatting and checks for element
#  definition existence and parameter definition existence by element. Does not do any conversions.
class WarpWriter:
    element_template = "{name} = wp.{type}("
    parameter_template = " {parameter}={value},"
    line_template = "{line_name} = "

    def __init__(self, beamline, input_simulation):
        self.beamline = beamline

        to_warp = 'to_madx'  # Just a reminder of what this is doing
        self.mapping = _build_field_map(_FIELD_MAP)[input_simulation][to_warp]

        self.element_names = []  # list of written elements
        self.drift_name = 'Drft'

    def write_file(self, filename):
        contents = "import warp as wp\n"
        line_contents = self.line_template.format(line_name=self.beamline.name)
        elements = self.beamline.get_beamline_elements()
        
        # Write element definitions
        for ele in elements:
            name = ele.name.upper()
            ele_type =  ele.type.upper()
            
            # Write element definition
            if ele_type in self.mapping.keys() and name not in self.element_names:
                new_ele_type = self.mapping[ele_type][0]
                ele_str = self.element_template.format(name=name, type=new_ele_type)
                for par in self.mapping[ele_type][1:]:
                    try:
                        ele_str += self.parameter_template.format(parameter=par, value=ele.parameters[par])
                    except KeyError:
                        raise KeyError(f'Could not find {par} in element {name} of type {ele_type}')
                contents += ele_str + ')' +'\n'
                
                self.element_names.append(name)
            elif name in self.element_names:
                continue
            else:
                print("Element '{name}' with type {type} is not defined. Will be added as drift.".format(name=name, type=ele.type))
                ele_str = self.element_template.format(name=name, type=self.drift_name)
                if ele.parameters.get('l'):
                    length = ele.parameters.get('l')
                else:
                    length = 0.0
                ele_str += self.parameter_template.format(parameter='l', value=length)

                contents += ele_str + ')' + '\n'
                
                self.element_names.append(name)
            # Append to line defition
            line_contents += '{} + '.format(name)
            
        contents += line_contents[:-2]  # Cut off trailing '+'
        
        with open(filename, 'w') as ff:
            ff.write(contents)
            
        return contents
