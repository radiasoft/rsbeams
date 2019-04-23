#====================================================================
#
# YuE comment (04/11/19):
#
# Script realize approach to simulate IOTA bunch dynamics with changing of 
# parameters of the nonlinear elements after each defined number of turns.
#
# File 'multuNLsimulation_example.py' from directory
#              /home/vagrant/jupyter/eidelyur/iota
# was used as source to create this script.
#
#====================================================================
#
# Import and setup IPython magics:
#
%matplotlib inline
%reload_ext autoreload
%autoreload 2
import sys, os
import numpy as np
import scipy
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec
#
# Synergia specific imports:
#
import rssynergia 
from rssynergia.base_diagnostics import read_bunch
from rssynergia.base_diagnostics import workflow
from rssynergia.base_diagnostics import lfplot
from rssynergia.base_diagnostics import plotbeam
from rssynergia.base_diagnostics import latticework
from rssynergia.base_diagnostics import basic_calcs
from rssynergia.base_diagnostics import pltbunch
from rssynergia.base_diagnostics import elliptic_sp
from rssynergia.base_diagnostics import options
from rssynergia.base_diagnostics import diagplot
from rssynergia.elliptic import elliptic_beam6d
import synergia
import synergia_workflow

def next_break():
    selection = 'loop'
    while selection == 'loop':
        selection = raw_input('\nYour selection: next, break')
        print 'You selectiion is ', selection
        if (selection == 'next'):
            return 1
        if (selection == 'break'):
            return 2

def plotcoordDistr(Particles,numbParticles):
    newCoordinates = np.zeros((6,numbParticles))
    for k in range(numbParticles):
        for j in range(6):
            newCoordinates[j,k] = Particles[0,k,j]        
    xmax = np.max(newCoordinates[0,:])
    xpmax = np.max(newCoordinates[1,:])
    ymax = np.max(newCoordinates[2,:])
    ypmax = np.max(newCoordinates[3,:])

    #another way - use gridspec
    fig = plt.figure(figsize=(15,5))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1,1,1]) 
    
    ax0 = plt.subplot(gs[0])
    plt.plot(newCoordinates[0,:],newCoordinates[1,:],'.',color='r')
    ax0.set_title('X-X\' Coordinate Space',fontsize='16')
    ax0.set_xlim([-1.5*xmax,1.5*xmax])
    ax0.set_ylim([-1.5*xpmax,1.5*xpmax])
    
    ax1 = plt.subplot(gs[1])
    plt.plot(newCoordinates[2,:],newCoordinates[3,:],'.',color='b')
    ax1.set_title('Y-Y\' Coordinate Space',fontsize='16')
    ax1.set_xlim([-1.5*ymax,1.5*ymax])
    ax1.set_ylim([-1.5*ypmax,1.5*ypmax])
    
    ax2 = plt.subplot(gs[2])
    plt.plot(newCoordinates[0,:],newCoordinates[2,:],'.',color='k')
    ax2.set_title('X-Y Trace Space',fontsize='16')
    ax2.set_xlim([-1.5*xmax,1.5*xmax])
    ax2.set_ylim([-1.5*ymax,1.5*ymax])
    
    fig.canvas.set_window_title('Synergia Phase Space Distribution')
    fig.tight_layout()
    plt.show()
    return
#
# Names of the structures:
#
nameStruct = []
numbStruct = 2
for n in range(numbStruct):
    nameStruct.append('nlnr'+str(n)+'_IO_82')
# Checking:
for n in range(numbStruct):
    print 'nameStrct[{}] = {}'.format(n,nameStruct[n])
#    
# Import of lattices from MADX files:
#
lattices = {}
# Dictionary of lattices:
lattice_repo = '/home/vagrant/jupyter/eidelyur/iota/'

for n in range(numbStruct):
    nameCrrnt = nameStruct[n]
    if (n == 0):
        lattices[nameCrrnt] = lattice_repo + "lattice_1IO_center.madx" #centered t1 8.2 1IO lattice
    if (n == 1):
#        lattices[nameCrrnt] = lattice_repo + "lattice_1IO_nll_center.madx" #centered t1 8.2 1IO lattice
        lattices[nameCrrnt] = lattice_repo + "lattice_1IO_center.madx" #centered t1 8.2 1IO lattice
# Checking:
for n in range(numbStruct):
    nameCrrnt = nameStruct[n]
    print 'lattices[{}] = {}'.format(nameCrrnt,lattices[nameCrrnt])
#
# Some global parameters:
#
nsteps_per_element = 5
order = 1
output_dir = 'example_run'
workflow.make_path(output_dir)
numbMacroPrtcls = 10000
tval = 0.4         # elliptic strength parameter
cval = 0.01        # aperture parameter
init_tune = 0.3
lnl_length = 1.8   # m
nSegmnt = 20
# We want the normalized emittance in x to be 0.3 mm-mrad
init_norm_emittance = 0.3*1.e-6
turns_per_structure = 100

#
# Some flags (plot, if =1):
#
plotBunchFlag = 1
plotCoordinatesFlag = 1
plotStructureFlag = 1
#
# Loop over all structures:
#
# for n in range(numbStruct):
for n in range(2):
    nameCrrnt = nameStruct[n]
    latticeCrrnt = synergia.lattice.MadX_reader().get_lattice("iota", lattices[nameCrrnt])

    for elem in latticeCrrnt.get_elements():
        if elem.get_type() == 'nllens':
            elem.set_string_attribute("extractor_type", "chef_propagate")
        else:
            elem.set_string_attribute("extractor_type", "chef_map") 
    nstepsCrrnt = len(latticeCrrnt.get_elements())*nsteps_per_element
    optsCrrnt = workflow.make_opts(nameCrrnt, order, output_dir, nstepsCrrnt, nsteps_per_element)
    optsCrrnt.macro_particles = numbMacroPrtcls
    optsCrrnt.steps_per_element = nsteps_per_element
    stepperCrrnt = synergia.simulation.Independent_stepper_elements \
                   (latticeCrrnt, optsCrrnt.map_order, optsCrrnt.steps_per_element)
    lattice_simulatorCrrnt = stepperCrrnt.get_lattice_simulator()
# Plotting of the structure:
    optsCrrnt.lattice_name = 'IOTA 8-2: file '+ lattices[nameCrrnt]
    optsCrrnt.ID = None
    optsCrrnt.path = None
    optsCrrnt.output_dir = output_dir
    optsCrrnt.turns = turns_per_structure 
    optsCrrnt.variance = 0.5
    optsCrrnt.lattice_simulator = lattice_simulatorCrrnt
    optsCrrnt.relpath = output_dir
    optsCrrnt.lf_fns = ['beta_x','beta_y','D_x']
    optsCrrnt.lattice = latticeCrrnt
    optsCrrnt.save = False
    optsCrrnt.scale = 2
    if (plotStructureFlag == 1):
        lfplot.plot_sliced_lattice_functions(optsCrrnt)
# Data to construct a matched bunch:    
    ref = latticeCrrnt.get_reference_particle()
    beta = ref.get_beta()
    gamma = ref.get_gamma()
    if n == 0:
        optsCrrnt.norm_emittance = init_norm_emittance
        optsCrrnt.emitx = basic_calcs.calc_geometric_emittance \
                          (optsCrrnt.norm_emittance, beta, gamma)
        optsCrrnt.emity = optsCrrnt.emitx
# Construct a matched bunch:
        myBunch = synergia.optics.generate_matched_bunch_transverse \
                  (lattice_simulatorCrrnt, optsCrrnt.emitx, optsCrrnt.emity, optsCrrnt.stdz, \
                   optsCrrnt.dpop, optsCrrnt.real_particles,optsCrrnt.macro_particles, optsCrrnt.seed) 
# Plot a matched bunch:
        if (plotBunchFlag == 1):
            pltbunch.plot_bunch(myBunch)      # Plotting of x-y, x-px and y-py distributions
# Data to calculate the parameters of nonlinear element::
    optsCrrnt.t = tval                # will use for any n
    optsCrrnt.c = cval                # will use for any n
    optsCrrnt.lnll = lnl_length       # will use for any n
    optsCrrnt.nseg = nSegmnt          # will use for any n
    optsCrrnt.new_tune = init_tune    # will use for n=0 only!
#        
# Construct the nonlinear element
# Attention: for n > 1 it is necessary to redefine the parameter optsCrrnt.new_tune;
# Method returns: f0, betae,alphae,betas:
# f0 - focal length, betae - entrance beta function, alphae - entrance alpha function,
# betas - middle beta function
#
    vals = basic_calcs.get_base_nll(optsCrrnt.lnll, optsCrrnt.new_tune, optsCrrnt.t, optsCrrnt.c)
#
# Using of parameters of nonlinear element to generate a toy costing beam with fixed 
# elliptical Hamiltonian and returns the particle array, while also saving the array to
# a text file. Coordinates are chosen with fixed random number generator seed, so that they
# should always produce the same initial distribution for a given emittance.
#   
    optsCrrnt.betae = vals[3]
    optsCrrnt.alphae = 0 #fixed 0 alpha for center
    optsCrrnt.beta0 = vals[3]
    optsCrrnt.emits = [9.74e-6]
    optsCrrnt.lattice = latticeCrrnt
    optsCrrnt.quiet = False  
    particlesCrrnt = elliptic_beam6d.toyellipticalbeam6D(optsCrrnt)
# Plot a toy beam (x-x', y-y', x-y coordinates):
    if (plotCoordinatesFlag == 1):
        plotcoordDistr(particlesCrrnt,numbMacroPrtcls)
# Construct a toyheader for quick calculation of bunch properties
    toyheader = {}
    toyheader['s_val'] = 0.
    optsCrrnt.lattice_simulator = lattice_simulatorCrrnt  
    for index in range(len(optsCrrnt.emits)):
        bunchCrrnt = particlesCrrnt[index]
        initialH,initialI = elliptic_sp.calc_bunch_H(bunchCrrnt,optsCrrnt)
        bunch_mean = np.mean(initialH)
        bunch_std = np.std(initialH)
        bunch_var = (bunch_std/bunch_mean)*100
        print "Constructed bunch with {} macroparticles, having mean H: {} and std: {}%". \
              format(optsCrrnt.macro_particles, bunch_mean,bunch_var)
# now add longitudinal momentum variation
# For random samples with mean = 0, sigma = sigma, use sigma*np.random.randn(...)
        bunchCrrnt[:,5] = optsCrrnt.dpop*np.random.randn(1,len(bunchCrrnt))
# bunch[:,5] = np.zeros(len(bunch)) #0 dpop

        optsCrrnt.num_total_particles = optsCrrnt.macro_particles*len(optsCrrnt.emits)
        optsCrrnt.tracked_particles = optsCrrnt.num_total_particles
    particles_file = '{}/myBunch.txt'.format(optsCrrnt.output_dir)

    np.savetxt(particles_file,bunchCrrnt)            #  write the bunch to a text file
    bucket_length = beta*latticeCrrnt.get_length()/4 # RF harmonic number is 4

    commCrrnt = synergia.utils.Commxx(True)          # define a communicator
    myBunch = read_bunch.read_bunch(particles_file,ref,optsCrrnt.real_particles,commCrrnt,bucket_length)

    if (plotBunchFlag == 1):
        pltbunch.plot_bunch(myBunch)      # Plotting of x-y, x-px and y-py distributions
#============================#
#                            #
# Perform a basic simulation #
#                            #
#============================#

# We will run our matched beam through the nonlinear lattice for 10 turns for each 
# structure, outputing individual particle coordinates (Diagnostics_particles) each turn
# and basic RMS bunch properties (Diagnostics_basic) each step (slice) of the simulation.

turns_per_structure = 100

#
# Construct a matched bunch (before first structure):
#
nameFirstStrctr = nameStruct[0]
latticeFirstStrctr = synergia.lattice.MadX_reader().get_lattice("iota", lattices[nameFirstStrctr])
for elem in latticeFirstStrctr.get_elements():
    if elem.get_type() == 'nllens':
        elem.set_string_attribute("extractor_type", "chef_propagate")
    else:
        elem.set_string_attribute("extractor_type", "chef_map") 
nstepsFirstStrctr = len(latticeFirstStrctr.get_elements())*nsteps_per_element
optsFirstStrctr = workflow.make_opts(nameFirstStrctr,order,output_dir,nstepsFirstStrctr,nsteps_per_element)
optsFirstStrctr.macro_particles = numbMacroPrtcls
optsFirstStrctr.steps_per_element = nsteps_per_element
stepperMain = synergia.simulation.Independent_stepper_elements \
              (latticeFirstStrctr, optsFirstStrctr.map_order, optsFirstStrctr.steps_per_element)
lattice_simulatorMain = stepperMain.get_lattice_simulator()
# Plotting of the structure:
optsFirstStrctr.lattice_name = 'IOTA 8-2: file '+ lattices[nameFirstStrctr]
optsFirstStrctr.ID = None
optsFirstStrctr.path = None
optsFirstStrctr.turns = turns_per_structure 
optsFirstStrctr.variance = 0.5
optsFirstStrctr.lattice_simulator = lattice_simulatorMain
optsFirstStrctr.relpath = output_dir
optsFirstStrctr.lf_fns = ['beta_x','beta_y','D_x']
optsFirstStrctr.lattice = latticeFirstStrctr
optsFirstStrctr.save = False
optsFirstStrctr.scale = 2
lfplot.plot_sliced_lattice_functions(optsFirstStrctr)
#
# Data to construct this bunch and buiding it:    
#
ref = latticeFirstStrctr.get_reference_particle()
beta = ref.get_beta()
gamma = ref.get_gamma()
optsFirstStrctr.norm_emittance = init_norm_emittance
optsFirstStrctr.emitx = basic_calcs.calc_geometric_emittance \
                        (optsFirstStrctr.norm_emittance, beta, gamma)
optsFirstStrctr.emity = optsFirstStrctr.emitx
comm = synergia.utils.Commxx(True) #define a communicator
myBunch = synergia.optics.generate_matched_bunch_transverse \
          (lattice_simulatorMain,optsFirstStrctr.emitx, optsFirstStrctr.emity, optsFirstStrctr.stdz, \
           optsFirstStrctr.dpop, optsFirstStrctr.real_particles,optsFirstStrctr.macro_particles, \
           optsFirstStrctr.seed) 
# Plot this bunch:
pltbunch.plot_bunch(myBunch)

bunch_simulator = synergia.simulation.Bunch_simulator(myBunch)

print "Before main simulation"
for n in range(1):
    nameCrrnt = nameStruct[n]
    latticeCrrnt = synergia.lattice.MadX_reader().get_lattice("iota", lattices[nameCrrnt])
    for elem in latticeCrrnt.get_elements():
        if elem.get_type() == 'nllens':
            elem.set_string_attribute("extractor_type", "chef_propagate")
        else:
            elem.set_string_attribute("extractor_type", "chef_map") 
    nstepsCrrnt = len(latticeCrrnt.get_elements())*nsteps_per_element
    optsCrrnt = workflow.make_opts(nameCrrnt, order, output_dir, nstepsCrrnt, nsteps_per_element)
    optsCrrnt.macro_particles = numbMacroPrtcls
    optsCrrnt.steps_per_element = nsteps_per_element
    print "Before stepperCrrnt"
    stepperCrrnt = synergia.simulation.Independent_stepper_elements \
                   (latticeCrrnt, optsCrrnt.map_order, optsCrrnt.steps_per_element)
    print "Before lattice_simulatorCrrnt"
    lattice_simulatorCrrnt = stepperCrrnt.get_lattice_simulator()
    print "After lattice_simulatorCrrnt"

# basic diagnostics - PER STEP
#    crrntRept = synergia.bunch.Diagnostics_basic.get_repetition(myBunch)
#    crrnt_s = synergia.bunch.Diagnostics_basic.get_s(myBunch)
#    crrnt_s_n = synergia.bunch.Diagnostics_basic.get_s_n(myBunch)
#    print ('n = ',n,': crrntRept = ',crrntRept)
#    print ('crrnt_s = ',crrnt_s)
#    print ('crrnt_s_n = ',crrnt_s_n)
    basicFileName = 'basic_strct'+str(int(n))+'.h5'
    basicdiag = synergia.bunch.Diagnostics_basic(basicFileName, optsCrrnt.output_dir)
    bunch_simulator.add_per_step(basicdiag)
    print "Saving basic diagnostics each step to file {}".format(basicFileName)
    nextBreakFlag = int(next_break())
    if (nextBreakFlag == 1):
        num_partcls = int(synergia.bunch.Diagnostics_basic.get_num_particles())
        print 'num_partcls =', num_partcls
    
# include full diagnostics
    full2FileName = 'full2_strct'+str(int(n))+'.h5'
    fulldiag = synergia.bunch.Diagnostics_full2(full2FileName, optsCrrnt.output_dir)
    bunch_simulator.add_per_turn(fulldiag)
    print "Saving full2 diagnostics each turn to file {}".format(full2FileName)

# particle diagnostics - PER TURN
    optsCrrnt.turnsPerDiag = 10
    particleFileName = 'particles_strct'+str(int(n))+'.h5'
    particlediag = synergia.bunch.Diagnostics_particles(particleFileName,0,0,optsCrrnt.output_dir)
    bunch_simulator.add_per_turn(particlediag, optsCrrnt.turnsPerDiag)
    print "Saving turn-by-turn particle data every {} turns to file {}". \
          format(optsCrrnt.turnsPerDiag,particleFileName)

# track diagnostics - PER TURN
#    trackFileName = 'tracks_strct'+str(int(n))+'.h5'
#    try:
#        trackdiag = synergia.bunch.Diagnostics_track(trackFileName,0,optsCrrnt.output_dir)
#        print 'No problems with trackdiag...'
#    except:
#        print 'Problems with trackdiag...'
#    optsCrrnt.turnsPerDiag = [1]
#    try:
#        bunch_simulator.add_per_turn(trackdiag, optsCrrnt.turnsPerDiag)
#        print "Saving turn-by-turn track data every {} turns to file {}". \
#              format(optsCrrnt.turnsPerDiag,trackFileName)
#    except:
#        print 'Problems with saving of tracks...'

    optsCrrnt.turns = 100
    optsCrrnt.checkpointperiod = 10
    optsCrrnt.maxturns = optsCrrnt.turns+1

    myrank = comm.get_rank()
    print "setting up propagator for rank {}".format(myrank)
    print "Before propagator"
    propagator = synergia.simulation.Propagator(stepperCrrnt)
    print "After propagator and before checkpoint_priod"
    propagator.set_checkpoint_period(optsCrrnt.checkpointperiod)

    print "Starting simulation for rank {}".format(myrank)
    propagator.propagate(bunch_simulator,optsCrrnt.turns, optsCrrnt.maxturns,optsCrrnt.verbosity)
    print "End of simulation"

# clean up files
    workflow.cleanup(optsCrrnt.output_dir)

  