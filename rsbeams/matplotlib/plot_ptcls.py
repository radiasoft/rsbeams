# -*- coding: utf-8 -*-
"""Generalized algorithm for plotting contour and/or scatter plots.

Adapted from open source method: scatter_contour.py
https://github.com/astroML/astroML/blob/master/astroML/plotting/scatter_contour.py

Parameters
----------
plot_flag : style of plot (scatter, contour, line, etc.)
plot_type : axis scaling (linear, log-log, or semi-log)
x_data    : x data for the contour plot
y_data    : y data for the contour plot
ax        : the axes on which to plot
divs      : desired number of divisions along each axis
levels    : integer or array (optional, default=10)
            number of contour levels, or array of contour levels

Original code taken from RadTrack project, https://github.com/radiasoft/radtrack
:copyright: Copyright (c) 2013 RadiaBeam Technologies, LLC. All Rights Reserved.

Subsequent mods are due to RadiaSoft,
:copyright: Copyright (c) 2017 Radiasoft LLC. All Rights Reserved.

:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""

import math
import numpy
from matplotlib.path import Path

def scatter_contour(plot_flag, plot_type, x_data, y_data, ax, divs=10, levels=10):
    """Generate the specified plot"""
    ref = None
    if plot_flag in ['contour', 'combo']:
        if type(x_data) is list: # x_data contains data for 2 axis ranges
            levels = numpy.asarray(levels)
            if levels.size == 1:
                levels = numpy.linspace(min(y_data), max(y_data), levels)

            min_x = min(x_data[0])
            max_x = max(x_data[0])

            min_y = min(x_data[1])
            max_y = max(x_data[1])

            points = len(y_data)
            ratio = float(max_x - min_x)/(max_y - min_y)
            shape_x = math.sqrt(points*ratio)
            shape_y = math.sqrt(points/ratio)
            x_mesh, y_mesh = numpy.meshgrid(x_data[0], x_data[1])
            z_mesh = y_data.reshape([shape_x, shape_y])
            z_mesh = z_mesh[0:len(x_data[0]), 0:len(x_data[1])]

            ref = ax.contourf(x_mesh, y_mesh, z_mesh, levels=levels, \
                              extent=[min_x, max_x, min_y, max_y])

        else:
            threshold = 8 if plot_flag == 'combo' else 1

            # generate the 2D histogram, allowing the algorithm to use
            #   all data points, automatically calculating the 2D extent
            my_hist, xbins, ybins = numpy.histogram2d(x_data, y_data, divs)

            # specify contour levels, allowing user to input simple integer
            levels = numpy.asarray(levels)
            # if user specified an integer, then populate levels reasonably
            if levels.size == 1:
                levels = numpy.linspace(threshold, my_hist.max(), levels)

            # define the 'extent' of the contoured area, using the
            #   the horizontal and vertical arrays generaed by histogram2d()
            extent = [xbins[0], xbins[-1], ybins[0], ybins[-1]]
            i_min = numpy.argmin(levels)

            # draw a zero-width line, which defines the outer polygon,
            #   in order to reduce the number of points drawn
            outline = ax.contour(my_hist.T, levels[i_min:i_min+1],linewidths=0,extent=extent)

            # generate the contoured image, filled or not
            #   use my_hist.T, rather than full my_hist, to limit extent of the contoured region
            #   i.e. only the high-density regions are contoured
            #   the return value is potentially useful to the calling method
            ref = ax.contourf(my_hist.T, levels, extent=extent)

    # logic for finding particles in low-density regions
    if plot_flag == 'combo':
        # create new 2D array that will hold a subset of the particles
        #   i.e. only those in the low-density regions
        low_density_arr = numpy.hstack([x_data[:, None], y_data[:, None]])

        # extract only those particles outside the high-density region
        if len(outline.allsegs[0]) > 0:
            outer_poly = outline.allsegs[0][0]
            points_inside = Path(outer_poly).contains_points(low_density_arr)
            x_plot = low_density_arr[~points_inside]
        else:
            x_plot = low_density_arr

    if plot_flag.startswith('scatter') or plot_flag.endswith('line'):
        x_plot = numpy.hstack([x_data[:, None], y_data[:, None]])

    if plot_flag in ['combo', 'scatter', 'scatter-line']:

        # Below is a hack to get around the fact that scatter plots
        # don't get correct axis limits if either axis is log scale.
        #
        # ax.plot(...) seems to work, so draw a plot and then delete
        # it, leaving the plot with a correct axes view.

        toRemove, = ax.plot(x_plot[:,0], x_plot[:,1], c='w')
        ax.scatter(x_plot[:,0], x_plot[:,1], marker='.', c='k')
        ax.lines.remove(toRemove)

    if plot_flag.endswith('line'):
        ax.plot(x_plot[:,0], x_plot[:,1], c='k')

    if plot_flag in ['line', 'scatter', 'scatter-line']:
        if plot_type in ['log-log', 'semi-logx']:
            ax.set_xscale('log', nonposx='mask')

        if plot_type in ['log-log', 'semi-logy']:
            ax.set_yscale('log', nonposy='mask')

        if plot_type in ['linear', 'semi-logy']:
            ax.set_xscale('linear')

        if plot_type in ['linear', 'semi-logx']:
            ax.set_yscale('linear')

    return ref


def gen_contour_levels(field, nLevels=40, multiplier=1.1):
    """generate the contour levels"""

    max_e = multiplier * numpy.max(field)
    min_e = multiplier * numpy.min(field)
    if abs(min_e) < max_e:
        max_e = numpy.around(max_e, decimals=3)
        min_e = -max_e
    else:
        min_e= numpy.around(min_e, decimals=3)
        max_e = abs(min_e)

    e_levels = []
    delta_e = (max_e-min_e) / nLevels
    for iLoop in range(nLevels):
        e_levels.append(min_e + iLoop*delta_e)

    return e_levels
