import numpy as np
import sys

#######################
# Physical Parameters #
#######################
Lx,Ly = (2*np.pi, 2*np.pi)

# Viscosity
nu = 2.0e-6
# Hyperviscosity exponent ( nu*nabla^(2 nn) v)
nn = 2

# Kinetic energy of initial condition:
E0 = 1e-1

# Forcing amplitude
fp0 = 0
ft0 = 0

# Reduction factor

r = 1.0

# Planetary Beta value beta<1 dry unstable, beta>1 dry stable

beta = 0.78

# boundary tilt, h = 0 untilted, h = 1 tilted

h = 0

# drag coefficient, drag = \sqrt(2)*R_{LH}/2

drag = 0.11

# Forcing wavenumbers
kf_min = 8
kf_max = 10

#########################

# Noies amplitude
NA = 0.1

########################
# Numerical Parameters #
########################
# Numerical resolution
Nx = 128
Ny = 128

# Set how long you want to run for:
sim_tmax = np.inf  # simulation time units
real_tmax = 3*60*60  #(11+50/60.)*60*60 # 12*60= 12 mins. Real time is in seconds ... 12 hours = 43200
# Simulation will stop when the first of either is reached.

# Real (wall) time interval between snapshots exporting
tsnap_wall = 10 #real_tmax/3.-15*60 # So that we save 3 snapshots 15 minutes before each third

# # Time steps between outputting spectra and fluxes:
# tspec = 230

# Time steps between outputting time series info
tseries = 1 #23
