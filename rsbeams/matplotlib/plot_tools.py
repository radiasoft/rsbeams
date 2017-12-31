import math
import numpy as np
from matplotlib.path import Path


"""
Generalized algorithm for plotting contour and/or scatter plots.
  plotFlag is queried to determine what's done.

Adapted from open source method: scatter_contour.py
https://github.com/astroML/astroML/blob/master/astroML/plotting/scatter_contour.py

Parameters
----------
plotFlag : style of plot (scatter, contour, line, etc.)
plotType : axis scaling (linear, log-log, or semi-log)
x, y     : x and y data for the contour plot
ax       : the axes on which to plot
divs     : desired number of divisions along each axis
levels   : integer or array (optional, default=10)
         number of contour levels, or array of contour levels
"""
def scatConPlot(plotFlag, plotType, x, y, ax, divs=10, levels=10):
    ref = None
    if plotFlag in ['contour', 'combo']:
        if type(x) is list: # x contains data for 2 axis ranges
            levels = np.asarray(levels)
            if levels.size == 1:
                levels = np.linspace(min(y), max(y), levels)

            minX = min(x[0])
            maxX = max(x[0])

            minY = min(x[1])
            maxY = max(x[1])

            points = len(y)
            ratio = float(maxX - minX)/(maxY - minY)
            shapeX = math.sqrt(points*ratio)
            shapeY = math.sqrt(points/ratio)
            X, Y = np.meshgrid(x[0], x[1])
            Z = y.reshape([shapeX, shapeY])
            Z = Z[0:len(x[0]), 0:len(x[1])]

            ref = ax.contourf(X, Y, Z, levels=levels, extent=[minX, maxX, minY, maxY])

        else:
            threshold = 8 if plotFlag == 'combo' else 1

            # generate the 2D histogram, allowing the algorithm to use
            #   all data points, automatically calculating the 2D extent
            myHist, xbins, ybins = np.histogram2d(x, y, divs)

            # specify contour levels, allowing user to input simple integer
            levels = np.asarray(levels)
            # if user specified an integer, then populate levels reasonably
            if levels.size == 1:
                levels = np.linspace(threshold, myHist.max(), levels)

            # define the 'extent' of the contoured area, using the
            #   the horizontal and vertical arrays generaed by histogram2d()
            extent = [xbins[0], xbins[-1], ybins[0], ybins[-1]]
            i_min = np.argmin(levels)

            # draw a zero-width line, which defines the outer polygon,
            #   in order to reduce the number of points drawn
            outline = ax.contour(myHist.T, levels[i_min:i_min+1],linewidths=0,extent=extent)

            # generate the contoured image, filled or not
            #   use myHist.T, rather than full myHist, to limit extent of the contoured region
            #   i.e. only the high-density regions are contoured
            #   the return value is potentially useful to the calling method
            ref = ax.contourf(myHist.T, levels, extent=extent)

    # logic for finding particles in low-density regions
    if plotFlag == 'combo':
        # create new 2D array that will hold a subset of the particles
        #   i.e. only those in the low-density regions
        lowDensityArray = np.hstack([x[:, None], y[:, None]])

        # extract only those particles outside the high-density region
        if len(outline.allsegs[0]) > 0:
            outer_poly = outline.allsegs[0][0]
            points_inside = Path(outer_poly).contains_points(lowDensityArray)
            Xplot = lowDensityArray[~points_inside]
        else:
            Xplot = lowDensityArray

    if plotFlag.startswith('scatter') or plotFlag.endswith('line'):
        Xplot = np.hstack([x[:, None], y[:, None]])

    if plotFlag in ['combo', 'scatter', 'scatter-line']:

        # Terrible hack to get around the "fact" that scatter plots
        # do not get correct axis limits if either axis is log scale.
        # ax.plot(...) seems to work, so draw a plot and then delete
        # it, leaving the plot with a correct axes view.

        toRemove, = ax.plot(Xplot[:,0], Xplot[:,1], c='w')
        ax.scatter(Xplot[:,0], Xplot[:,1], marker='.', c='k')
        ax.lines.remove(toRemove)

    if plotFlag.endswith('line'):
        ax.plot(Xplot[:,0], Xplot[:,1], c='k')

    if plotFlag in ['line', 'scatter', 'scatter-line']:
        if plotType in ['log-log', 'semi-logx']:
            ax.set_xscale('log', nonposx='mask')

        if plotType in ['log-log', 'semi-logy']:
            ax.set_yscale('log', nonposy='mask')

        if plotType in ['linear', 'semi-logy']:
            ax.set_xscale('linear')

        if plotType in ['linear', 'semi-logx']:
            ax.set_yscale('linear')

    return ref



# function to generate contour levels
def generateContourLevels(field, nLevels=40, multiplier=1.1):
    # generate symmetric min/max values
    eMax = multiplier * np.max(field)
    eMin = multiplier * np.min(field)
    if abs(eMin) < eMax:
        eMax = np.around(eMax, decimals=3)
        eMin = -eMax
    else:
        eMin= np.around(eMin, decimals=3)
        eMax = abs(eMin)
        
    # create the level values
    eLevels = []
    deltaE = (eMax-eMin) / nLevels
    for iLoop in range(nLevels):
        eLevels.append(eMin + iLoop*deltaE)

    return eLevels
