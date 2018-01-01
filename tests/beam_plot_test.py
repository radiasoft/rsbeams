from __future__ import absolute_import, division, print_function, unicode_literals
import pytest
import numpy
from rsbeams.physics import rsconst
from rsbeams.statistics import stats6d
from rsbeams.ptcl_beam import RsPtclBeam6D
from rsbeams.ptcl_beam import RsPhaseSpace6D
from rsbeams.matplotlib import RsPlotPhaseSpace6D

def test_beam_plot():

    # Instantiate an arbitrary e- bunch
    num_ptcls = 1000
    design_p_ev = 271e+6
    total_charge_c = 3.05e-09
    mass_ev = rsconst.m_e_EV

    dist_type = 'gaussian'
    max_rms_fac = 4.9

    alpha_x = 1.3     # []
    beta_x = 21.1     # [m/rad]
    emit_x = 1.7e-06  # [m-rad]

    alpha_y = -3.01   # []
    beta_y = 9.07     # [m/rad]
    emit_y = 0.44e-6  # [m-rad]

    alpha_z = 0       # []
    beta_z =  1.55    # [m/rad]
    emit_z = 3.1e-05  # [m-rad]

    my_ebeam = RsPtclBeam6D.RsPtclBeam6D(num_ptcls, design_p_ev, \
                                         total_charge_c, mass_ev, \
                                         dist_type, max_rms_fac, \
                                         alpha_x, beta_x, emit_x, \
                                         alpha_y, beta_y, emit_y, \
                                         alpha_z, beta_z, emit_z )

    # Grab a pointer to the phase space object (particle data)
    my_space = my_ebeam.get_distrib6d().get_phase_space_6d()

    # Instantiate a bunch plotting object
    my_plotter = RsPlotPhaseSpace6D.RsPlotPhaseSpace6D(my_space)

    # Generate an x-px scatter plot
    my_plotter.set_title('x-px projection')
    my_plotter.plot_data6d(0,1)

    # Generate an x-px scatter plot (after the drift)
    my_plotter.set_title('y-py projection')
    my_plotter.plot_data6d(2,3)

    # Generate a y-pz scatter plot (after the drift)
    my_plotter.set_title('y-pz projection')
    my_plotter.plot_data6d(2,5)

    # Render the plots
    my_plotter.show_plots()


test_beam_plot()
