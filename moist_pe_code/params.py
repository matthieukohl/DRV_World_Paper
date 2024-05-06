import numpy as np
import sys

#######################
# Physical Parameters #
#######################
Lx,Ly = (6*np.pi, 6*np.pi)
Lz = 1
# Viscosity
nu = 5.0e-5
nu_theta = 5.0e-5
# Hyperviscosity exponent ( nu*nabla^(2 nn) v)
nn = 2

# Kinetic energy of initial condition:
E0 = 1e-1

# Forcing amplitude
fp0 = 0
ft0 = 0

# Reduction factor

r = 0.01

# Planetary Beta value beta<1 dry unstable, beta>1 dry stable

beta = 0.78

# boundary tilt, h = 0 untilted, h = 1 tilted

h = 0

# epsilon (Rossby number)

eps = 0.4

# drag coefficient, drag = \sqrt(2)*R_{LH}/2

drag = 0.11

# relax coefficient

relax = 0.35

# large scale inverse diffusion

hnu = 5
hnn = -1;

# Forcing wavenumbers
kf_min = 8
kf_max = 10

#########################

# Noise amplitude
NA = 1.0

########################
# Numerical Parameters #
########################
# Numerical resolution
Nx = 128
Ny = 128
Nz = 10
# Set how long you want to run for:
sim_tmax = 160  # simulation time units
real_tmax = np.inf  #11*60*60  #(11+50/60.)*60*60 # 12*60= 12 mins. Real time is in seconds ... 12 hours = 43200
# Simulation will stop when the first of either is reached.

# Real (wall) time interval between snapshots exporting
#tsnap_wall = 1*60 #real_tmax/3.-15*60 # So that we save 3 snapshots 15 minutes before each third

# Sim time

tsnap_sim = 0.5

# # Time steps between outputting spectra and fluxes:
tspec = 230

# Time steps between outputting time series info
tseries = 23
