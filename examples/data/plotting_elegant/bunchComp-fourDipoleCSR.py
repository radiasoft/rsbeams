#!/usr/bin/env python

# bunchComp - fourDipoleCSR
execution_mode = 'parallel'

lattice_file = """

"Q": CHARGE,total=1e-9
"B1": CSRCSBEND,angle=0.146607657167524,bins=500,csr=0,e2=0.146607657167524,l=0.200718260855179,n_slices=100,output_file="B1.output_file.sdds",output_interval=10,sg_halfwidth=1
"B2": CSRCSBEND,angle=-0.146607657167524,bins=500,csr=0,e1=-0.146607657167524,l=0.200718260855179,n_slices=100,output_file="B2.output_file.sdds",output_interval=10,sg_halfwidth=1
"B3": CSRCSBEND,angle=-0.146607657167524,bins=500,csr=0,e2=-0.146607657167524,l=0.200718260855179,n_slices=100,output_file="B3.output_file.sdds",output_interval=4,sg_halfwidth=1
"B4": CSRCSBEND,angle=0.146607657167524,bins=500,csr=0,e1=0.146607657167524,l=0.200718260855179,n_slices=100,output_file="B4.output_file.sdds",output_interval=4,sg_halfwidth=1
"L1": CSRDRIFT,dz=0.01,l=0.758132998376353,use_stupakov=1
"L2": CSRDRIFT,dz=0.01,l=0.5,use_stupakov=1
"L3": CSRDRIFT,dz=0.01,l=1,use_stupakov=1
"PF": PFILTER,deltalimit=0.005
"LINA10": RFCA,change_p0=1,freq=2856e6,l=0.3,phase=62,volt=4800000
"LINB10": RFCA,change_p0=1,freq=2856e6,l=0.3,phase=96,volt=4800000
"ZWAKE": WAKE,factor=8.6,inputfile="WAKE-inputfile.knsl45.liwake.sdds",interpolate=1,n_bins=1024,tcolumn="t",wcolumn="W"
"W1": WATCH,filename="w1.filename.sdds"
"W2": WATCH,filename="w2.filename.sdds"
"W3": WATCH,filename="w3.filename.sdds"
"W4": WATCH,filename="w4.filename.sdds"
"W5": WATCH,filename="w5.filename.sdds"
"LINACA": LINE=("LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE","LINA10","ZWAKE")
"LINACB": LINE=("LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE")
"CHICANE": LINE=("Q","LINACA","W1","B1","L1","W2","B2","L2","W3","B3","L1","W4","B4","W5","L3")
"BL": LINE=("CHICANE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","LINB10","ZWAKE","PF")

"""

elegant_file = """

&global_settings
  mpi_io_write_buffer_size = 1048576,
&end

&run_setup
  semaphore_file = run_setup.semaphore,
  centroid = "run_setup.centroid.sdds",
  lattice = "elegant.lte",
  output = "run_setup.output.sdds",
  p_central_mev = 54.9835,
  parameters = "run_setup.parameters.sdds",
  print_statistics = 1,
  sigma = "run_setup.sigma.sdds",
  use_beamline = "bl",
&end

&run_control
&end

&twiss_output
  alpha_x = 1,
  alpha_y = 1,
  beta_x = 10,
  beta_y = 10,
  filename = "twiss_output.filename.sdds",
  matched = 0,
  output_at_each_step = 1,
  statistics = 1,
&end

&bunched_beam
  alpha_x = 1,
  alpha_y = 1,
  beta_x = 10,
  beta_y = 10,
  distribution_cutoff[0] = 3, 3, 3,
  emit_x = 4.6e-08,
  emit_y = 4.6e-08,
  enforce_rms_values[0] = 1, 1, 1,
  n_particles_per_bunch = 50000,
  one_random_bunch = 0,
  sigma_dp = 0.001,
  sigma_s = 0.00065,
  symmetrize = 1,
&end

&track
&end

"""

with open('elegant.lte', 'w') as f:
    f.write(lattice_file)

with open('elegant.ele', 'w') as f:
    f.write(elegant_file)

import os
os.system('elegant elegant.ele')
