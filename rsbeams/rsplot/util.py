from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import numpy as np

_element_mappings = {'QUAD': ['KQUAD', 'QUAD'],
                     'SBEND': ['SBEND', 'CSBEND', 'CSRCSBEND', 'RBEND'],
                     'SEXT': ['SEXT', ],
                     'RFCW': ['RFCW', 'RFCA']}
_element_colors = {'QUAD': 'tab:green',
                   'SBEND': 'tab:red',
                   'SEXT': 'tab:blue',
                   'RFCW': 'tab:orange'}


def _invert_mapping(mapping):
    # Inverts element_mappings
    # Expects dict val's to be lists

    inverted_mapping = {}
    for key, val in mapping.items():
        for el in val:
            inverted_mapping[el] = key

    return inverted_mapping


def plot_profile(sdds_columns, axes, scale=1.0, height=0.25):
    """
    Plot a beamline magnet profile using matplotlib.
    Args:
        sdds_columns: Structured NumPy array from an elegant file (easiest to load with rsbeams.rsdata.SDDS.readSDDS).
                      Should contained columns: 's' and 'ElementType'
        axes: matplotlib axes object that profile will be plotted on.
        scale: Physical length scale to use. Default of 1.0 -> meters.
        height: Height of magnet profile Rectangles. Can be used to adjust visual representation. (default=0.25)

    Returns: matplotlib.collections.PatchCollection of Rectangles representation magnets.

    """
    zp = 0.0  # Vertical midpoint
    used_element_types = _invert_mapping(_element_mappings)
    s = 0.0
    patch_collection = []
    colors = []

    # Place horizontal line under elements
    axes.plot([0.0, sdds_columns['s'][-1] * scale], [0.0, 0.0], c='k')

    for i in range(sdds_columns.shape[-1]):
        element_type = sdds_columns['ElementType'][i]

        if element_type in used_element_types.keys():
            base_type = used_element_types[element_type]
            e1, e2 = s * scale, sdds_columns['s'][i] * scale
            axes.plot([s, e1], [zp, zp], c=_element_colors[base_type])
            mag = Rectangle((e1, -height / 2.),
                            width=e2 - e1, height=height)
            patch_collection.append(mag)
            colors.append(_element_colors[base_type])
            s = sdds_columns['s'][i]
        else:
            s = sdds_columns['s'][i]

    collection = PatchCollection(patch_collection, facecolor=colors, edgecolor=colors)
    axes.add_collection(collection)

    return collection

def get_histogram_points(data, bin_number):
    counts, bins = np.histogram(data, bins=int(bin_number))
    central_bins = np.array([(bins[i + 1] + bins[i]) / 2. for i in range(len(bins) - 1)])

    return central_bins, counts

def get_current_from_t(time_data, charge):
    q_per_p = charge / time_data.size

    bin_number = 256

    counts, bins = np.histogram((time_data - np.average(time_data)), bins=int(bin_number))
    bin_length = bins[1] - bins[0]
    # print(np.std(time_data))
    central_bins = np.array([(bins[i + 1] + bins[i]) / 2. for i in range(len(bins) - 1)])
    normalized_counts = counts * q_per_p / bin_length  # [count] * [Coulomb / count] / [seconds]

    return central_bins, normalized_counts
