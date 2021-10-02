import numpy as np
import matplotlib.pyplot as plt
from rsbeams.rsplot import util


def beamline_profile(sdds, page=0, quantities=None, xlim=None, ylim=None, save=None):
    sdds_columns = sdds.columns[page]

    if not quantities:
        quantities = ['Sx', 'Sy']

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(18, 6),
                                   sharex=True,
                                   gridspec_kw={'height_ratios': [1, 8]})

    ax1.axis('off')
    util.plot_profile(sdds_columns, ax1, height=0.25)
    for quant in quantities:
        ax2.plot(sdds_columns['s'], sdds_columns[quant], label=util.format_symbol(sdds.column_symbol(quant)))

    ax2.legend(fontsize=16)
    ax2.set_xlabel('s (m)', fontsize=16)

    ylabel = ','.join([util.format_symbol(sdds.column_symbol(quant)) for quant in quantities])
    ylabel += ' (' + ','.join([util.format_symbol(sdds.column_units(quant)) for quant in quantities]) + ')'
    ax2.set_ylabel(ylabel, fontsize=16)
    if xlim:
        ax2.set_xlim(*xlim)
    if ylim:
        ax2.set_ylim(*ylim)
    if save:
        plt.savefig(save)
    plt.show()


def phase_space(sdds, page=0, bins=128, save=None):
    """
    Create density plot of common projections: (x, x'), (y, y'), (x, y), (t, delta)
    Default units are mm and mrad for distances and angles; ps for time; and fractional deviation (delta) for momentum.
    Args:
        sdds: Open SDDS file from rsbeams.rsdata.readSDDS
        page: [0] Optional page number to read from
        bins: [128] Integer value for bins to create histogram
        save: [None] Full file name, including extension, if image should be saved

    Returns:

    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    labels = [('x (mm)', 'x\' (mrad)'),
              ('y (mm)', 'y\' (mrad)'),
              ('x (mm)', 'y (mm)'),
              ('t (ps)', r'$\delta$ ( )')]

    sdds_columns = sdds.columns[page]
    sdds_parameters = sdds.parameters[page]

    data = [(sdds_columns['x'] * 1e3, sdds_columns['xp'] * 1e3),
            (sdds_columns['y'] * 1e3, sdds_columns['yp'] * 1e3),
            (sdds_columns['x'] * 1e3, sdds_columns['y'] * 1e3),
            ((sdds_columns['t'] - np.average(sdds_columns['t'])) * 1e12,
             (sdds_columns['p'] - np.average(sdds_columns['p'])) / np.average(sdds_columns['p']))]

    for ax, datum, label in zip(axes.flatten(), data, labels):
        counts, xbins, ybins = np.histogram2d(*datum, bins=bins)
        ax.imshow(counts.T, interpolation='gaussian',
                  extent=(np.min(xbins), np.max(xbins), np.min(ybins), np.max(ybins)),
                  origin='lower',
                  aspect='auto')
        ax.set_xlabel(label[0])
        ax.set_ylabel(label[1])

    if save:
        plt.savefig(save)

    plt.show()


def _get_phase_space_figure(xs=12, ys=12):
    fig = plt.figure(figsize=(xs, ys))
    gs = fig.add_gridspec(3, 3)
    hm = fig.add_subplot(gs[:2, :2])
    mom = fig.add_subplot(gs[:2, -1])
    tim = fig.add_subplot(gs[-1, :2])

    return fig, gs, hm, mom, tim


def longitudinal_phase_space(sdds, page=0, charge='auto', save=None):
    """
    Plot of the longitudinal phase space with projection plots of current distribution and momentum profile.
    Args:
        sdds: Open SDDS file from rsbeams.rsdata.readSDDS
        page: [0] Optional page number to read from
        charge: ('auto') if auto the charge will be read from the SDDS file. If the charge cannot be found current
         profile will just be normalized to peak intensity at 1. If a float value is given this will be used to set the
         current (should be given in units of Coulombs)
        save: [None] Full file name, including extension, if image should be saved

    Returns:

    """
    fig, gs, hm, mom, tim = _get_phase_space_figure()

    sdds_columns = sdds.columns[page]
    sdds_parameters = sdds.parameters[page]

    time_scale = 1e12
    t = (sdds_columns['t'] - np.average(sdds_columns['t']))
    p = sdds_columns['p']

    normalized_charge = False
    if charge == 'auto':
        try:
            charge = sdds_parameters['Charge']
            if charge == 0.0:
                print('Warning Charge is 0 in SDDS file. Using normalized counts.')
                charge = 1.0
                normalized_charge = True
        except KeyError:
            print('Warning: No charge data found in SDDS file. Using normalized counts.')
            print('Set charge manually if desired.')
            charge = 1.0
            normalized_charge = True

    time_bins, current_counts = util.get_current_from_t(t, charge)
    momentum_bins, momentum_counts = util.get_histogram_points(p, 256)
    if normalized_charge:
        current_counts = current_counts / np.max(current_counts)
        time_label = 'Counts ()'
    else:
        time_label = 'Current (A)'

    counts, xbins, ybins = np.histogram2d(t * time_scale, p, bins=128)
    #     im = ax.imshow(H.T, cmap=cmap, interpolation='sinc', norm=matplotlib.colors.LogNorm())

    hm.imshow(counts.T, interpolation='gaussian',
              extent=(np.min(xbins), np.max(xbins), np.min(ybins), np.max(ybins)),
              origin='lower',
              aspect='auto')
    hm.set_ylabel('p ($m_e c$)')
    hm.set_xticklabels([])

    mom.plot(momentum_counts, momentum_bins)
    mom.set_yticklabels([])
    mom.set_xlabel('Counts ()')

    tim.plot(time_bins * time_scale, current_counts)
    tim.set_xlabel('t (ps)')
    tim.set_ylabel(time_label)

    if save:
        plt.savefig(save)

    plt.show()


def horizontal_phase_space(sdds, page=0, save=None):
    phase_space_projection(sdds, 'x', 'xp', x_label='x (mm)', y_label='x\' (mrad)', page=page, save=save)


def vertical_phase_space(sdds, page=0, save=None):
    phase_space_projection(sdds, 'y', 'yp', page=page, x_label='y (mm)', y_label='y\' (mrad)', save=save)


def phase_space_projection(sdds, coordinate_x, coordinate_y, x_label='', y_label='', page=0, save=None):
    """
    Plot of the longitudinal phase space with projection plots of current distribution and momentum profile.
    Args:
        sdds: Open SDDS file from rsbeams.rsdata.readSDDS
        page: [0] Optional page number to read from
        save: [None] Full file name, including extension, if image should be saved

    Returns:

    """
    fig, gs, hm, mom, tim = _get_phase_space_figure()

    sdds_columns = sdds.columns[page]

    space_scale = 1e3
    x = sdds_columns[coordinate_x]
    y = sdds_columns[coordinate_y]


    x_bins, x_counts = util.get_histogram_points(x * space_scale, 256)
    y_bins, y_counts = util.get_histogram_points(y * space_scale, 256)
    counts, mxbins, mybins = np.histogram2d(x * space_scale, y * space_scale, bins=128)

    hm.imshow(counts.T, interpolation='gaussian',
              extent=(np.min(mxbins), np.max(mxbins), np.min(mybins), np.max(mybins)),
              origin='lower',
              aspect='auto')

    hm.set_ylabel(y_label)
    hm.set_xticklabels([])

    mom.plot(y_counts, y_bins)
    mom.set_yticklabels([])
    mom.set_xlabel('Counts ()')

    tim.plot(x_bins, x_counts)
    tim.set_xlabel(x_label)
    tim.set_ylabel('Counts ()')

    if save:
        plt.savefig(save)

    plt.show()