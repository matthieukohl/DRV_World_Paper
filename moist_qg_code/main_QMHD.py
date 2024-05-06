"""
Dedalus script for 3D, rotating quasistatic magnetohydrodynamics (QMHD). 
This regime is valid in the limit of low magnetic Reynolds number and 
a strong uniform, background mangetic field.
This script uses a Fourier basis in all directions with periodic boundary
conditions.  
"""

import numpy as np
import h5py
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

import time
import pathlib
from dedalus import public as de
from dedalus.extras import flow_tools

import logging
logger = logging.getLogger(__name__)

ncores = comm.Get_size()
nnodes = int(ncores/16)
logger.info('Running on %s cores,  %s nodes' % (ncores,nnodes))

from params import *

from fractional_laplacian import FractionalLaplacian

#############################
# DEFINITIONS AND OPERATORS #
#############################
import spec1d
import flux1d
import nonlin_term
import vector_cal as vc

# Add operator to namespace
de.operators.parseables['FractionalLaplacian'] = FractionalLaplacian

##########
# DOMAIN #
##########

# Create bases and domain
x_basis = de.Fourier('x', Nx, interval=(0, Lx), dealias=3/2)
y_basis = de.Fourier('y', Ny, interval=(0, Ly), dealias=3/2)
z_basis = de.Fourier('z', Nz, interval=(0, Lz), dealias=3/2)

domain = de.Domain([x_basis, y_basis, z_basis], grid_dtype=np.float64,mesh=(nnodes,16))

# For general use
x = domain.grid(0)
y = domain.grid(1)
z = domain.grid(2)
kx = domain.elements(0)
ky = domain.elements(1)
kz = domain.elements(2)
k2 = kx**2 + ky**2 + kz**2

###########
# FORCING #
###########
# Initialize fields
Fx = domain.new_field()
Fy = domain.new_field()
Fz = domain.new_field()

# Forcing range
cond = (k2<=kf_max**2)&(k2>=kf_min**2)

local_coeff_shape=Fx['c'].shape
phase = np.random.uniform(low=-np.pi,high=np.pi,size=local_coeff_shape)
Fx['c'][cond] = (np.cos(phase[cond])+1j*np.sin(phase[cond]))/np.sqrt(k2[cond])

local_coeff_shape=Fy['c'].shape
phase = np.random.uniform(low=-np.pi,high=np.pi,size=local_coeff_shape)
Fy['c'][cond] = (np.cos(phase[cond])+1j*np.sin(phase[cond]))/np.sqrt(k2[cond])

local_coeff_shape=Fz['c'].shape
phase = np.random.uniform(low=-np.pi,high=np.pi,size=local_coeff_shape)
Fz['c'][cond] = (np.cos(phase[cond])+1j*np.sin(phase[cond]))/np.sqrt(k2[cond])

# Make sure it's divergence-free
Fx2,Fy2,Fz2 = vc.curl(domain,[Fx,Fy,Fz])
Fx['c']=Fx2['c']
Fy['c']=Fy2['c']
Fz['c']=Fz2['c']

# Normalize:
KE_op = 0.5*(Fx*Fx + Fy*Fy + Fz*Fz)/(Lx*Ly*Lz)
KE_int_op = de.operators.integrate(KE_op,'x','y','z') 
E = KE_int_op.evaluate()['g'][0,0,0]
Fx['g'] = Fx['g']*f0/np.sqrt(E)
Fy['g'] = Fy['g']*f0/np.sqrt(E)
Fz['g'] = Fz['g']*f0/np.sqrt(E)

#############
# EQUATIONS #
#############
# 3D QMHD
problem = de.IVP(domain, variables=['p','u','v','w'])
# Rotation/Magnetic field
problem.parameters['NNx'] = NNx
problem.parameters['NNz'] = NNz
problem.parameters['Omega'] = Omega
# Dissipation
problem.parameters['nu'] = nu
problem.substitutions['nn'] = "%f" % nn # Hyperviscosity exponent
problem.parameters['hnu'] = hnu
problem.substitutions['mm'] = "%f" % mm # Hypoviscosity exponent
# Forcing
problem.parameters['Fx'] = Fx
problem.parameters['Fy'] = Fy
problem.parameters['Fz'] = Fz
# Operator substitutions
problem.substitutions['Adv(a)'] = "u*dx(a) + v*dy(a) + w*dz(a)"
problem.substitutions['FL'] = "FractionalLaplacian"
problem.substitutions['BO(a)'] = "FL(NNx**2*dx(dx(a)) + 2*NNx*NNz*dx(dz(a)) + NNz**2*dz(dz(a)),-1)"
# Equtions of motion
problem.add_equation("dx(u) + dy(v) + dz(w) = 0", condition="(nx != 0) or (ny != 0) or (nz != 0)")
problem.add_equation("p = 0", condition="(nx == 0) and (ny == 0) and (nz == 0)")
                              #Viscosity           #Rotation  #JxB    #Hypoviscoity
problem.add_equation("dt(u) + nu*FL(u,nn) + dx(p) - Omega*v + BO(u) + hnu*FL(u,-mm)  = -Adv(u) + Fx")
problem.add_equation("dt(v) + nu*FL(v,nn) + dy(p) + Omega*u + BO(v) + hnu*FL(v,-mm)  = -Adv(v) + Fy")
problem.add_equation("dt(w) + nu*FL(w,nn) + dz(p)           + BO(w) + hnu*FL(w,-mm)  = -Adv(w) + Fz")

# Build solver
solver = problem.build_solver(de.timesteppers.RK222)
#de.timesteppers.SBDF3)
#de.timesteppers.RK443)

logger.info('Solver built')

#################################
# Initial conditions or restart #
#################################
if not pathlib.Path('restart.h5').exists():

    # Initial conditions
    u = solver.state['u']
    v = solver.state['v']
    w = solver.state['w']
    p = solver.state['p']

    # Random initial conditions for k <= kf_max, otherwise = 0
    cond = k2>kf_max**2

    local_coeff_shape=u['c'].shape
    phase = np.random.uniform(low=-np.pi,high=np.pi,size=local_coeff_shape)
    u['c'][k2!=0] = (np.cos(phase[k2!=0])+1j*np.sin(phase[k2!=0]))/np.sqrt(k2[k2!=0])
    u['c'][cond] = 0.0
    u['c'][k2==0] = 0.0

    local_coeff_shape=v['c'].shape
    phase = np.random.uniform(low=-np.pi,high=np.pi,size=local_coeff_shape)
    v['c'][k2!=0] = (np.cos(phase[k2!=0])+1j*np.sin(phase[k2!=0]))/np.sqrt(k2[k2!=0])
    v['c'][cond] = 0.0
    v['c'][k2==0] = 0.0

    local_coeff_shape=w['c'].shape
    phase = np.random.uniform(low=-np.pi,high=np.pi,size=local_coeff_shape)
    w['c'][k2!=0] = (np.cos(phase[k2!=0])+1j*np.sin(phase[k2!=0]))/np.sqrt(k2[k2!=0])
    w['c'][cond] = 0.0
    w['c'][k2==0] = 0.0

    # Make sure it's divergence-free
    u2,v2,w2 = vc.curl(domain,[u,v,w])
    u['c']=u2['c']
    v['c']=v2['c']
    w['c']=w2['c']

    # Normalize:
    KE_op = 0.5*(u*u + v*v + w*w)/(Lx*Ly*Lz)
    KE_int_op = de.operators.integrate(KE_op,'x','y','z') 
    E = KE_int_op.evaluate()['g'][0,0,0]
    u['g'] = u['g']*np.sqrt(E0/E)
    v['g'] = v['g']*np.sqrt(E0/E)
    w['g'] = w['g']*np.sqrt(E0/E)

    # Calculate the initial pressure based on these velocity fields
    p2 = nonlin_term.pressure_calc(domain,[u,v,w])
    p['c']=p2['c']

    # Timestepping and output
    dt = 1e-5 # Make initial guess (extra small, just in case)
    fh_mode = 'overwrite'
    
else:
    # Restart
    write, last_dt = solver.load_state('restart.h5', -1)

    # Timestepping and output
    dt = last_dt
    fh_mode = 'append'

# Integration parameters
solver.stop_sim_time = sim_tmax
solver.stop_wall_time = real_tmax

#######
# CFL #
#######
CFL = flow_tools.CFL(solver, initial_dt=dt, safety=0.5, cadence=10, threshold=0.05)
CFL.add_velocities(('u','v','w'))
# Rotation rate
if Omega>0:
        CFL.add_frequency(Omega)
        logger.info('1/dt restriction from rotation = %f' % Omega)

# Momentum diffusion
kcut = Nx/2. # not Nx/3 because of the 3/2 dealiasing rule instead of the 2/3 rule.
CFL.add_frequency(nu*kcut**(2*nn))
CFL.add_frequency(hnu)

# B_0 operator frequency
CFL.add_frequency(NNx**2+2*NNx*NNz+NNz**2)
    
logger.info('1/dt restriction from visc = %f, hypovisc = %f, and B0 = %f' % (nu*kcut**(2*nn),hnu,NNx**2+2*NNx*NNz+NNz**2))

############
# Analysis #
############
# SNAPSHOTS
snapshots = solver.evaluator.add_file_handler('snapshots', wall_dt=tsnap_wall, max_writes=50, mode=fh_mode)
snapshots.add_system(solver.state)

# TIME SERIES
t_series = solver.evaluator.add_file_handler('time_series', iter=tseries,mode=fh_mode)
t_series.add_task("integ(u*u+v*v+w*w)/2", name='en')
t_series.add_task("integ(u*Fx+v*Fy+w*Fz)",name = 'inj')
t_series.add_task("integ(u*nu*FL(u,nn)+v*nu*FL(v,nn)+w*nu*FL(w,nn))",name = 'diss')
t_series.add_task("integ(u*hnu*FL(u,-mm)+v*hnu*FL(v,-mm)+w*hnu*FL(w,-mm))",name = 'hdiss')
t_series.add_task("integ(u*BO(u) + v*BO(v) + w*BO(w))", name='bdiss')

# Flow properties # HELP: What is this exactly? What's the difference between flow properties and the tasks above?
flow = flow_tools.GlobalFlowProperty(solver, cadence=10)
flow.add_property("sqrt(u*u + v*v + w*w)", name='Re')
flow.add_property("(u*u+v*v+w*w)/2", name='KE')

# Initialized h5 datasets for flux and spectra (unless they already exist):
if rank==0:
    if fh_mode=='overwrite':
        f = h5py.File('flux_spec.h5','w')
        spec = f.create_group('spectra')
        flux = f.create_group('flux')
    else:
        f = h5py.File('flux_spec.h5','a')
        spec = f['spectra']
        flux = f['flux']
        
# Main loop
try:
    logger.info('Starting loop')
    start_time = time.time()
    while solver.proceed:
        dt = CFL.compute_dt()
        dt = solver.step(dt)
        if (solver.iteration-1) % 10 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
            logger.info('Max KE = %f' %flow.max('Re'))
            logger.info('KE = %f' %flow.volume_average('KE'))
        # Spectra and fluxes
        if (solver.iteration-1) % tspec == 0: 
            logger.info('Calculating fluxes and spectra,Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
            # Load fields
            u = solver.state['u']
            v = solver.state['v']
            w = solver.state['w']
            p = solver.state['p']
            
            # Calculate KE spectrum
            [modn,En] = spec1d.spec1d(domain,[u,v,w])
            
            # Calculate nonlinear term:
            a,b,c = nonlin_term.nonlin_term(domain,[u,v,w,p])
            
            # Calculate KE flux:
            [modn_f,Pi] = flux1d.flux1d_classic(domain,[u,v,w],[a,b,c],flux=True)
        
            # Save data:
            if rank==0:
                if ((solver.iteration-1)==0)&(fh_mode=='overwrite'):
                    spec.create_dataset('time', data = [solver.sim_time],maxshape=(None,),chunks=True)
                    spec.create_dataset('modn',data=modn)
                    spec.create_dataset('KE',data=[En],maxshape=(None,np.shape(En)[0]),chunks=True)
                    flux.create_dataset('time', data = [solver.sim_time],maxshape=(None,),chunks=True)
                    flux.create_dataset('modn',data=modn_f)
                    flux.create_dataset('KE',data=[Pi],maxshape=(None,np.shape(Pi)[0]),chunks=True)
                else:
                    spec['time'].resize((spec['time'].shape[0] + 1), axis = 0)
                    spec['time'][-1:] = [solver.sim_time]
                    spec['KE'].resize((spec['KE'].shape[0] + 1), axis = 0)
                    spec['KE'][-1:] = [En]
                    flux['time'].resize((flux['time'].shape[0] + 1), axis = 0)
                    flux['time'][-1:] = [solver.sim_time]
                    flux['KE'].resize((flux['KE'].shape[0] + 1), axis = 0)
                    flux['KE'][-1:] = [Pi]
        
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    if rank==0:
        f.close()
    end_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_time-start_time))
    logger.info('Run time: %f cpu-hr' %((end_time-start_time)/60/60*domain.dist.comm_cart.size))
