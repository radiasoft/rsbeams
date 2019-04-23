import os
import synergia
import numpy as np
import h5py
#import tables
#import argparse
import inspect
from mpi4py import MPI

# load the particles that will be used for the simulation
# The particles file is a text file with particle coordinates
# defined with the MAD-X conventions: X PX Y PY T PT
# Read this in using numpy's loadtxt command
# particle coordinates are converted to Synergia conventions

# input arguments:
#    particles_file: the file name
#    reference particle: the lattice reference particle for kinematic conversions
#    num_real_particles: the number of real particles in the bunch
#    bucket_length: the longitudinal length of the bucket
#    comm: the Commxx communicator object for this bunch
#    verbose: be chatty about what's happening
#  
def read_bunch(particles, refpart, real_particles, comm, bucket_length = None, verbose=False):
    '''
    Read a bunch from file (either .txt, .h5, or .mxtxt (MAD-X txt file)) and construct a Synergia bunch object.
    
    Arguments:
        - particles (string or np.ndarray):                        EITHER a file containing particles coordinates OR an ndarray of coordinates
        - refpart (synergia.foundation.foundation.Reference_particle):  the Synergia reference particle describing the bunch
        - num_real_particles (float):                                   the number of real particles
        - comm (synergia.utils.parallel_utils.Commxx):                  the Commxx communicator object for this bunch
        - bucket_length (Optional[float]):                              if specified, the longitudinal length of the bucket in m
        - verbose (Optional[Boolean]):                                  Flag for verbose feedback
    
    Returns:
        -bunch: A Synergia bunch object is created in the current session
    '''

    #first attempt to load the particles as an h5 file       
    try:
        return read_h5_particles(particles, refpart, real_particles, bucket_length, comm, verbose)
    
    #it's not an h5 file - then there are two possibilities:
    #1. It's another sort of file, in which case, an IOError will be thrown
    #2. It's a numpy array, in which case a TypeError will be thrown
    #Therefore, we will catch the IOErrror and process it as an input file to check if it's a legible text file
    #Then we will catch the possible TypeError and process it for being a numpy array
    
    except IOError:
        #IOError, so it's a file but not an .h5 file
        name,extension = os.path.splitext(particles)
            
        #assuming no error is thrown, we continue processing the file - whihc should be now either a .txt or .mxtxt
        assert extension == '.txt' or extension == '.mxtxt', \
        "Supported file types are hdf5 (.h5) and plain text (.txt/.mxtxt)"
            
        return read_txt_particles(particles, refpart, real_particles, bucket_length, comm, extension == '.mxtxt', verbose)
                    
    except TypeError:
        #TypeError, so it's not a file - so we should check if it's a numpy array
        #Had we checked the .txt read first, it would have return an AttributeError
        assert isinstance(particles, np.ndarray), \
        "Supported data types are numpy arrays only."
            
        return read_array_particles(particles, refpart, real_particles, bucket_length, comm, verbose)

#====================================================================

# if madx_format is True, the particles are in madX units, otherwise they are in
# synergia units

def read_txt_particles(particles_file, refpart, real_particles, bucket_length, comm, madx_format, verbose):
    """Read an array of particles from a text file"""
    
    four_momentum = refpart.get_four_momentum()
    pmass = four_momentum.get_mass()
    E_0 = four_momentum.get_total_energy()
    p0c = four_momentum.get_momentum()

    myrank = comm.get_rank()
    mpisize = comm.get_size()
    
    if myrank==0 and verbose:
        if madx_format:
            print "Loading madX particles from txt file: ", particles_file
        else:
            print "Loading Synergia particles from txt file: ", particles_file

    if myrank == 0:
        particles = np.loadtxt(particles_file)
        num_total_particles = particles.shape[0]
        # broadcast num particles to all nodes
        MPI.COMM_WORLD.bcast(num_total_particles, root=0)
    else:
        num_total_particles = None
        num_total_particles = MPI.COMM_WORLD.bcast(num_total_particles, root=0)

    if myrank == 0:
        # make sure the data has the correct shape, either [n,6] without
        # particles IDs or [n,7] with particle IDs.
        if (particles.shape[1] != 6) and (particles.shape[1] != 7):
            raise RuntimeError, "input data shape %shas incorrect number of particle coordinates"%repr(particles.shape)
        
        
        if madx_format:
            # numpy manipulations to convert kinematics
            # convert MAD-X T=-c*dt to Synergia c*ct
            particles[:,4] = -particles[:,4]
            # convert MAD-X Delta-E/pc to Synergia delta-p/p
            # sqrt(((dE/p0c)+(E0/p0c))**2 - (m/p0c)**2) - (p0c/p0c)
            m_over_pc = pmass/p0c
            E_0_over_pc = E_0/p0c
            particles[:,5] = np.sqrt( (particles[:,5] + E_0_over_pc) *
                                      (particles[:,5] + E_0_over_pc) - m_over_pc**2 ) - 1.0
        

        # if there are no IDs, append particle ID column
        if particles.shape[1] != 7:
            particles_w_id = np.column_stack((particles,
                                              np.arange(num_total_particles, dtype='d')))
        else:
            particles_w_id = particles
            
            if myrank == 0:
                print "Read ", num_total_particles, " particles"
    
    #Note: Synergia bunch constructor updated - commit 077b99d7 - 11/17/2016
    #Using old constructor throws an ArgumentError of a non-standard type.
    # Using a try and except to handle both instances.
    try:
        # try the original constructor
        bunch = synergia.bunch.Bunch(
            refpart,
            num_total_particles, real_particles, comm,
            bucket_length)
    except Exception, e:
        #look to see if it's an ArgumentError by evaluating the traceback
        if (not str(e).startswith("Python argument types in")):
            raise
        else:
            # use the new constructor
            if verbose:
                print "Using updated bunch constructor"
            bunch = synergia.bunch.Bunch(
                refpart,
                num_total_particles, real_particles, comm)
            # now set the new parameter 'z_period_length'
            if bucket_length is not None:
                bunch.set_z_period_length(bucket_length)
            else:
                bucket_length = 1. #fix this quantity

    local_num = bunch.get_local_num()
    local_particles = bunch.get_local_particles()

    # Each processor will have a possibly different number of local particles.
    # rank 0 has to find out how many each of them has and distribute them
    n_particles_by_proc = MPI.COMM_WORLD.gather(local_num, 0)
    if myrank == 0:
        # copy in my particles
        this_rank_start = 0
        local_particles[:,:] = particles_w_id[0:local_num, :]
        this_rank_start += local_num
        # send particles out to other ranks
        for r in range(1, mpisize):
            this_rank_end = this_rank_start+n_particles_by_proc[r]
            MPI.COMM_WORLD.send(obj=particles_w_id[this_rank_start:this_rank_end, :],
                                dest=r)
            this_rank_start += n_particles_by_proc[r]
    else:
        # I'm not rank 0.  Receive my particles
        lp = MPI.COMM_WORLD.recv(source=0)
        local_particles[:,:] = lp[:,:]
    return bunch

#==========================================================

def read_h5_particles(particles_file, refpart, real_particles, bucket_length, comm, verbose):
    """Read an array of particles from an HDF-5 file"""
    
    four_momentum = refpart.get_four_momentum()
    pmass = four_momentum.get_mass()
    E_0 = four_momentum.get_total_energy()
    p0c = four_momentum.get_momentum()

    myrank = comm.get_rank()
    mpisize = comm.get_size()
    
    if myrank==0 and verbose:
        print "Loading particles from h5 file: ", particles_file

    if myrank == 0:
        #h5 = tables.open_file(particles_file)
        h5 = h5py.File(particles_file)
        
        # use explicit int conversion otherwise there seems to
        # be a typepython->C++ type  mismatch of numpy.int64->int
        #num_total_particles = int(h5.root.particles.shape[0])
        num_total_particles = int(h5['particles'].shape[0])
        
        if verbose:
            print "Total of  ", num_total_particles, " particles from file"
        # broadcast num particles to all nodes
        MPI.COMM_WORLD.bcast(num_total_particles, root=0)
    else:
        num_total_particles = None
        num_total_particles = MPI.COMM_WORLD.bcast(num_total_particles, root=0)

    if myrank == 0:
        particles = h5['particles']
        # make sure the data has the correct shape, either [n,6] without
        # particles IDs or [n,7] with particle IDs.
        if (particles.shape[1] != 7):
            raise RuntimeError, "input data shape %shas incorrect number of particle coordinates"%repr(particles.shape)
    
    #Note: Synergia bunch constructor updated - commit 077b99d7 - 11/17/2016
    #Using old constructor throws an ArgumentError of a non-standard type.
    # Using a try and except to handle both instances.
    try:
        # try the original constructor
        bunch = synergia.bunch.Bunch(
            refpart,
            num_total_particles, real_particles, comm,
            bucket_length)
    except Exception, e:
        #look to see if it's an ArgumentError by evaluating the traceback
        if (not str(e).startswith("Python argument types in")):
            raise
        else:
            # use the new constructor
            if verbose:
                print "Using updated bunch constructor"
            bunch = synergia.bunch.Bunch(
                refpart,
                num_total_particles, real_particles, comm)        
            # now set the new parameter 'z_period_length'
            if bucket_length is not None:
                bunch.set_z_period_length(bucket_length)
            else:
                bucket_length = 1. #fix this quantity
            

    local_num = bunch.get_local_num()
    local_particles = bunch.get_local_particles()

    # Each processor will have a possibly different number of local particles.
    # rank 0 has to find out how many each of them has and distribute them
    n_particles_by_proc = MPI.COMM_WORLD.gather(local_num, 0)
    if myrank == 0:
        # copy in my particles
        this_rank_start = 0
        local_particles[:,:] = particles[0:local_num, :]
        this_rank_start += local_num
        # send particles out to other ranks
        for r in range(1, mpisize):
            this_rank_end = this_rank_start+n_particles_by_proc[r]
            MPI.COMM_WORLD.send(obj=particles[this_rank_start:this_rank_end, :],
                                dest=r)
            this_rank_start += n_particles_by_proc[r]
    else:
        # I'm not rank 0.  Receive my particles
        lp = MPI.COMM_WORLD.recv(source=0)
        local_particles[:,:] = lp[:,:]

    return bunch


#==========================================================

def read_array_particles(particle_array, refpart, real_particles, bucket_length, comm, verbose):
    """Read an array of particles coordinates from memory"""
    
    four_momentum = refpart.get_four_momentum()
    pmass = four_momentum.get_mass()
    E_0 = four_momentum.get_total_energy()
    p0c = four_momentum.get_momentum()

    myrank = comm.get_rank()
    mpisize = comm.get_size()
    
    if myrank==0 and verbose:
        print "Loading particles from: ".format(particle_array)

    if myrank == 0:
        
        # use explicit int conversion otherwise there seems to
        # be a typepython->C++ type  mismatch of numpy.int64->int
        #num_total_particles = int(h5.root.particles.shape[0])
        num_total_particles = particle_array.shape[0]
        
        if verbose:
            print "Total of  ", num_total_particles, " particles"
        # broadcast num particles to all nodes
        MPI.COMM_WORLD.bcast(num_total_particles, root=0)
    else:
        num_total_particles = None
        num_total_particles = MPI.COMM_WORLD.bcast(num_total_particles, root=0)

    if myrank == 0:
        particles = particle_array
        # make sure the data has the correct shape, either [n,6] without
        # particles IDs or [n,7] with particle IDs.
        if (particle_array.shape[1] != 7):
            raise RuntimeError, "input data shape %shas incorrect number of particle coordinates"%repr(particles.shape)
    
    #Note: Synergia bunch constructor updated - commit 077b99d7 - 11/17/2016
    #Using old constructor throws an ArgumentError of a non-standard type.
    # Using a try and except to handle both instances.
    try:
        # try the original constructor
        bunch = synergia.bunch.Bunch(
            refpart,
            num_total_particles, real_particles, comm,
            bucket_length)
    except Exception, e:
        #look to see if it's an ArgumentError by evaluating the traceback
        if (not str(e).startswith("Python argument types in")):
            raise
        else:
            # use the new constructor
            if verbose:
                print "Using updated bunch constructor"
            bunch = synergia.bunch.Bunch(
                refpart,
                num_total_particles, real_particles, comm)        
            # now set the new parameter 'z_period_length'
            if bucket_length is not None:
                bunch.set_z_period_length(bucket_length)
            else:
                bucket_length = 1. #fix this quantity
            

    local_num = bunch.get_local_num()
    local_particles = bunch.get_local_particles()

    # Each processor will have a possibly different number of local particles.
    # rank 0 has to find out how many each of them has and distribute them
    n_particles_by_proc = MPI.COMM_WORLD.gather(local_num, 0)
    if myrank == 0:
        # copy in my particles
        this_rank_start = 0
        local_particles[:,:] = particle_array[0:local_num, :]
        this_rank_start += local_num
        # send particles out to other ranks
        for r in range(1, mpisize):
            this_rank_end = this_rank_start+n_particles_by_proc[r]
            MPI.COMM_WORLD.send(obj=particles[this_rank_start:this_rank_end, :],
                                dest=r)
            this_rank_start += n_particles_by_proc[r]
    else:
        # I'm not rank 0.  Receive my particles
        lp = MPI.COMM_WORLD.recv(source=0)
        local_particles[:,:] = lp[:,:]

    return bunch



#================================================================

def print_bunch_stats(bunch):
    coord_names = ("x", "xp", "y", "yp", "c*dt", "dp/p")

    myrank = bunch.get_comm().get_rank()
    means = synergia.bunch.Core_diagnostics().calculate_mean(bunch)
    stds = synergia.bunch.Core_diagnostics().calculate_std(bunch, means)
    if myrank == 0:
        print "%20s   %20s   %20s"%("coord","mean","rms")
        print "%20s   %20s   %20s"%("====================",
                                    "====================",
                                    "====================")
        for i in range(6):
            print "%20s   %20.12e   %20.12e"%(coord_names[i], means[i], stds[i])

#=========================================================
