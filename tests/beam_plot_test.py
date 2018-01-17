from __future__ import absolute_import, division, print_function, unicode_literals

import os
import numpy
from rsbeams.physics import rsconst
from rsbeams.statistics import stats6d
from rsbeams.ptcl_beam import RsPtclBeam6D
from rsbeams.ptcl_beam import RsPhaseSpace6D
from rsbeams.matplotlib import RsScatterPS6D

def save_plot_and_check(plot_obj, file_name, force_error=False):
    """Encapsulate some file handling to avoid code repetition"""

    # File cleanup
    try:
        os.remove(file_name)
    except OSError:
        pass

    # save the plot to specified file name
    plot_obj.save_plots(file_name)

    # clear plotting object for next use
    plot_obj.clear_plots()

    # delete file for no good reason (testing the test)
    if (force_error == True):
        try:
            os.remove(file_name)
        except OSError:
            pass

    # is plot there with expected size?
    if os.path.exists(file_name):
#        print('File "',file_name,'" does exist.')
#        print(os.stat(file_name))
#        print('size: ', os.stat(file_name).st_size)
        assert(os.stat(file_name).st_size > 29000)
        assert(os.stat(file_name).st_size < 45000)
    else:
        # File doesn't exist; something went wrong
        assert 0, '\n   File "'+file_name+'" does NOT exist.'

    return


def test_beam_plot():
    """Test the plotting of beams"""

    # Decide if you want interactive plotting
    # Always use 'is_interactive = False' for testing
    is_interactive = False

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
    my_plotter = RsScatterPS6D.RsScatterPS6D(my_space)

    # Generate an x-px scatter plot
    my_plotter.set_title('x-xp projection')
    my_plotter.plot_data6d(0,1)
    if (is_interactive == False):
        save_plot_and_check(my_plotter, 'beam_plot_test_xxp.png')

    # Generate an x-px scatter plot (after the drift)
    my_plotter.set_title('y-yp projection')
    my_plotter.plot_data6d(2,3)
    if (is_interactive == False):
        save_plot_and_check(my_plotter, 'beam_plot_test_yyp.png')

    # Generate a y-pz scatter plot (after the drift)
    my_plotter.set_title('y-dp/p projection')
    my_plotter.plot_data6d(2,5)
    if (is_interactive == False):
#        save_plot_and_check(my_plotter, 'beam_plot_test_ydpop.png', True)
        save_plot_and_check(my_plotter, 'beam_plot_test_ydpop.png')

    # Render all plots to the screen
    if (is_interactive == True):
        my_plotter.show_plots()

# test_beam_plot()
