#!/usr/bin/env python

import synergia_workflow

opts = synergia_workflow.Options("iota")


opts.add("map_order", 1, "Map order", int)
#default output directory
opts.add("output_dir","SC_test","Directory for output files", str)
#opts.add("map_order", 1, "Map order", int)
#opts.add("steps", steps, "Number of steps per turn", int)
opts.add("steps_per_element",2,"Number of steps per element", int)


opts.add("verbosity", 1, "Verbosity of propagation", int)
opts.add("turns", 200, "Number of turns", int)
opts.add("maxturns", 2000, "Maximum number of turns to run before checkpointing and quitting", int)
opts.add("checkpointperiod", 3000, "Number of turns to run between checkpoints", int)


#opts.add("emitx", 2.5e-6, "real sigma Horizontal emittance [m rad]", float)
#opts.add("emity", 2.5e-6, "real sigma Vertical emittance [m rad]", float)
#opts.add("emit_transverse", 7.0e-6, "transverse emittance for elliptical beam [m rad]", float)
opts.add("radius", 0.5, "aperture radius [m]", float)
opts.add("emit",9.74e-6, "H0 value corresponding to real sigma horizontal emittance of 0.3 mm-mrad", float)
opts.add("stdz", 0.05, "sigma read z [m]", float) #5 cm bunch length for IOTA
opts.add("dpop", 0.0, "Delta-p/p spread", float)

opts.add("macro_particles", 10000, "Number of macro particles", int)
opts.add("real_particles", 1.0e11, "Number of real particles", float)
opts.add("tracked_particles", 10000, "Number of tracked particles", int)
opts.add("seed", 349250524, "Pseudorandom number generator seed", int)

opts.add("bunch_file","myBunch.txt","txt file for bunch particles", str)

#----------Space Charge Stuff---------------------
opts.add("gridx", 64, "grid points in x for solver", int)
opts.add("gridy", 64, "grid points in y for solver", int)
opts.add("gridz", 1, "grid points in z for solver", int) #1 for explicit 2D solver
opts.add("spacecharge", True, "whether space charge is on", bool)
opts.add("solver", "2dopen-hockney", "solver to use, '2dopen-hockney','3dopen-hockney', '2dbassetti-erskine'", str)

#options for controlling chef propagation vs. chef mapping!
opts.add("use_maps", "none", "use maps for propagation either all, none, onlyrf, nonrf")    #none means chef propagate
#opts.add("use_maps", "all", "use maps for propagation either all, none, onlyrf, nonrf")
#opts.add("allmaps", False, "Use all maps for propagation", bool)
opts.add("requested_stepper", "splitoperator", "Simulation stepper, either 'independent','elements','splitoperator','soelements'", str)
# YuE correction (02/20/19):
opts.add("stepper", "splitoperator", "Simulation stepper, either 'independent','elements','splitoperator','soelements'", str)
# opts.add("steps", steps, "Number of steps per turn", int)
# End of YuE correction

#----------MPI STUFF---------------------
opts.add("comm_divide", 8, "size of communicator")
#opts.add("matching", "4dmoments", "matching procedure (4dmoments)")
#opts.add("checkpointperiod", 2000, "Number of turns to run between checkpoints", int)
opts.add("concurrent_io", 8, "number of concurrent io threads for checkpointing", int)
#opts.add("nsigma", 8.0, "nsigma for solver", float)
#opts.add("long_kicks", True, "use longitudinal kicks", bool)
#opts.add("cutoffnsigma", 2.6, "cut off bunch at cutoffnsigma standard deviations")

#job_mgr = synergia_workflow.Job_manager("sc_test_14mA.py", opts)

#job_mgr = synergia_workflow.Job_manager("iota_66_1IO_nll_space_charge.py", opts)
