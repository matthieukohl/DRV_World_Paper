"""
Dedalus script for 2-layer moist QG equation.
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

domain = de.Domain([x_basis, y_basis], grid_dtype=np.float64,mesh=(nnodes,16))
# Domain with nnodes commented out

#domain = de.Domain([x_basis, y_basis], grid_dtype=np.float64)

# For general use
x = domain.grid(0)
y = domain.grid(1)
X,Y = np.meshgrid(x,y)
kx = domain.elements(0)
ky = domain.elements(1)
k2 = kx**2 + ky**2

###########
# FORCING #
###########
# Initialize fields
Fp = domain.new_field() # Phi frocing
Ft = domain.new_field() # Tau forcing

# Forcing range
cond = (k2<=kf_max**2)&(k2>=kf_min**2)

local_coeff_shape=Fp['c'].shape
phase = np.random.uniform(low=-np.pi,high=np.pi,size=local_coeff_shape)
Fp['c'][cond] = (np.cos(phase[cond])+1j*np.sin(phase[cond]))/np.sqrt(k2[cond])  # Dividing by |k| because the curl will give velocity

local_coeff_shape=Ft['c'].shape
phase = np.random.uniform(low=-np.pi,high=np.pi,size=local_coeff_shape)
Ft['c'][cond] = (np.cos(phase[cond])+1j*np.sin(phase[cond]))/np.sqrt(k2[cond])

# Normalize:
KE_op = (Fp*Fp)/(Lx*Ly) # Calculating the volume average of |Fp|^2
KE_int_op = de.operators.integrate(KE_op,'x','y')
E = KE_int_op.evaluate()['g'][0,0]
Fp['g'] = Fp['g']*fp0/np.sqrt(E)

KE_op = (Ft*Ft)/(Lx*Ly) # Calculating the volume average of |Fp|^2
KE_int_op = de.operators.integrate(KE_op,'x','y')
E = KE_int_op.evaluate()['g'][0,0]
Ft['g'] = Ft['g']*ft0/np.sqrt(E)


##########################
# Other fields/functions #
##########################

# Height field
#h = domain.new_field()
#h['g'] = 1 #Y

# Definition of r
def thresh_r(*args):
    #x   = args[0].data # this is an array; we use .data to get its values
    #y   = args[1].data
    w = args[0].data
    rtemp = args[1].value

    R = np.copy(w)
    R[w>=0]= rtemp
    R[w<0] = 1

    return R

def Sthresh_r(*args, domain=domain, F=thresh_r):
    """This function takes arguments *args, a function F, and a domain and
    returns a Dedalus GeneralFunction that can be applied.

    """
    return de.operators.GeneralFunction(domain, layout='g', func=thresh_r, args=args)

# Make it parseable
de.operators.parseables['thresh_r'] = Sthresh_r

#############
# EQUATIONS #
#############
problem = de.IVP(domain, variables=['phi','tau','w'])
# reduction factor
problem.parameters['r'] = r
# boundary tilt
problem.parameters['h'] = h
# beta value
problem.parameters['beta'] = beta
# drag value
problem.parameters['drag'] = drag
# Dissipation
problem.parameters['nu'] = nu
problem.substitutions['nn'] = "%f" % nn # Hyperviscosity exponent
# Forcing
problem.parameters['Fp'] = Fp
problem.parameters['Ft'] = Ft
# Other fields/functions
problem.parameters['h'] = h
# Operator substitutions
### Maybe we could define r here?
problem.substitutions['J(a,b)'] = "dx(a)*dy(b)-dy(a)*dx(b)"
problem.substitutions['Lap(a)'] = "dx(dx(a)) + dy(dy(a))"
problem.substitutions['Adv(a)'] = "dx(dx(dx(a))+dy(dy(a)))"
problem.substitutions['FL'] = "FractionalLaplacian"
# Equtions of motion
#                     cfl           +  nu*k**2*nn              k^2*phi            k^2*tau^2/phi    tau h/phi
problem.add_equation("dt(Lap(phi)) + nu*FL(Lap(phi),nn) + Adv(tau) + beta*dx(phi) - h*dx(tau) + drag*(Lap(phi)-Lap(tau)) =  -J(phi,Lap(phi)) - J(tau,Lap(tau))", condition="(nx != 0) or (ny != 0)")
problem.add_equation("dt(Lap(tau)) + nu*FL(Lap(tau),nn) + Adv(phi) + beta*dx(tau) - h*dx(phi) - drag*(Lap(phi)-Lap(tau)) + w = -J(phi,Lap(tau)) - J(tau,Lap(phi))", condition="(nx != 0) or (ny != 0)")
problem.add_equation(" Lap(w) - w + drag*(Lap(phi)-Lap(tau)) - 2*Adv(phi) - beta*dx(tau) + h*dx(phi) = Lap((1-thresh_r(w,r))*w) + 2*J(tau,Lap(phi)) - 2*J(dx(phi),dx(tau)) - 2*J(dy(phi),dy(tau))", condition="(nx != 0) or (ny != 0)") 
problem.add_equation("phi = 0" , condition="(nx == 0) and (ny == 0)")
problem.add_equation("tau = 0" , condition="(nx == 0) and (ny == 0)")
problem.add_equation("w = 0" , condition="(nx == 0) and (ny == 0)")

# Build solver
solver = problem.build_solver(de.timesteppers.RK443)
#de.timesteppers.SBDF3)
#de.timesteppers.RK443)
#de.timesteppers.RK222)

logger.info('Solver built')

# Build solver
solver = problem.build_solver(de.timesteppers.RK443)
#de.timesteppers.SBDF3)
#de.timesteppers.RK443)
#de.timesteppers.RK222)

logger.info('Solver built')

#################################
# Initial conditions or restart #
#################################
if not pathlib.Path('restart.h5').exists():

    # Initial conditions
    phi = solver.state['phi']
    tau = solver.state['tau']
    w = solver.state['w']

    local_coeff_shape=w['g'].shape
    tau['g'] = NA*np.random.uniform(size=local_coeff_shape)
    phi['g'] = NA*np.random.uniform(size=local_coeff_shape)
    w['g'] = NA*np.random.uniform(size=local_coeff_shape)

    # Take away small-scale noise
    tau['c'][k2>3**2] = 0.0
    w['c'][k2>3**2] = 0.0
    phi['c'][k2>3**2] = 0.0

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
CFL.add_velocities(('dx(phi)','dy(phi)'))
CFL.add_velocities(('dx(tau)','dy(tau)'))
#CFL.add_velocity('w')

# Momentum diffusion
kcut = Nx/2. # not Nx/3 because of the 3/2 dealiasing rule instead of the 2/3 rule.
CFL.add_frequency(nu*kcut**(2*nn))

## ADD other nonlinear terms

logger.info('1/dt restriction from visc = %f' % (nu*kcut**(2*nn)))

############
# Analysis #
############
# SNAPSHOTS
snapshots = solver.evaluator.add_file_handler('snapshots', wall_dt=tsnap_wall, max_writes=50, mode=fh_mode)
snapshots.add_system(solver.state)
snapshots.add_task("Lap(phi)" ,name='zeta_bt')
snapshots.add_task("Lap(tau)" ,name='zeta_bc')

## TIME SERIES
t_series = solver.evaluator.add_file_handler('time_series', iter=tseries,mode=fh_mode)
t_series.add_task("integ(dx(phi)**2 + dy(phi)**2 + dx(tau)**2 + dy(tau)**2)", name='en')
t_series.add_task("integ(w**2)", name='en_w')
#t_series.add_task("integ(u*Fx+v*Fy+w*Fz)",name = 'inj')
#t_series.add_task("integ(u*nu*FL(u,nn)+v*nu*FL(v,nn)+w*nu*FL(w,nn))",name = 'diss')
#t_series.add_task("integ(u*BO(u) + v*BO(v) + w*BO(w))", name='bdiss')
t_series.add_task("integ(phi*nu*FL(Lap(phi),nn)+phi*drag*(Lap(phi)-Lap(tau)))",name = 'phidiss')
t_series.add_task("integ(tau*nu*FL(Lap(tau),nn)-tau*drag*(Lap(phi)-Lap(tau)))",name = 'taudiss')

# Flow properties # HELP: What is this exactly? What's the difference between flow properties and the tasks above?
flow = flow_tools.GlobalFlowProperty(solver, cadence=10)
flow.add_property("sqrt(dx(phi)**2 + dy(phi)**2 + dx(tau)**2 + dy(tau)**2 +w*w)", name='Re')
flow.add_property("(dx(phi)**2 + dy(phi)**2 + dx(tau)**2 + dy(tau)**2 +w*w)/2", name='KE')

## Initialized h5 datasets for flux and spectra (unless they already exist):
if rank==0:
    if fh_mode=='overwrite':
        f = h5py.File('flux_spec.h5','w')
        spec = f.create_group('spectra')
        #flux = f.create_group('flux')
    else:
        f = h5py.File('flux_spec.h5','a')
        spec = f['spectra']
        #flux = f['flux']

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
            phi  = solver.state['phi']
            tau  = solver.state['tau']
            w = solver.state['w']
            phix = x_basis.Differentiate(phi)
            phiy = y_basis.Differentiate(phi)
            taux = x_basis.Differentiate(tau)
            tauy = y_basis.Differentiate(tau)

            # Calculate KE spectrum
            [modn,En_phi] = spec1d.spec1d(domain,[phix,phiy])
            [modn,En_tau] = spec1d.spec1d(domain,[taux,tauy])
            [modn,En_w] = spec1d.spec1d(domain,[w])
            # Calculate nonlinear term:
            #a,b,c = nonlin_term.nonlin_term(domain,[u,v,w,p])

            # Calculate KE flux:
            #[modn_f,Pi] = flux1d.flux1d_classic(domain,[u,v,w],[a,b,c],flux=True)

            # Save data:
            if rank==0:
                if ((solver.iteration-1)==0)&(fh_mode=='overwrite'):
                    spec.create_dataset('time', data = [solver.sim_time],maxshape=(None,),chunks=True)
                    spec.create_dataset('modn',data=modn)
                    spec.create_dataset('KE_phi',data=[En_phi],maxshape=(None,np.shape(En_phi)[0]),chunks=True)
                    spec.create_dataset('KE_tau',data=[En_tau],maxshape=(None,np.shape(En_tau)[0]),chunks=True)
                    spec.create_dataset('KE_w',data=[En_w],maxshape=(None,np.shape(En_w)[0]),chunks=True)
                    #flux.create_dataset('time', data = [solver.sim_time],maxshape=(None,),chunks=True)
                    #flux.create_dataset('modn',data=modn_f)
                    #flux.create_dataset('KE',data=[Pi],maxshape=(None,np.shape(Pi)[0]),chunks=True)
                else:
                    spec['time'].resize((spec['time'].shape[0] + 1), axis = 0)
                    spec['time'][-1:] = [solver.sim_time]
                    spec['KE_phi'].resize((spec['KE_phi'].shape[0] + 1), axis = 0)
                    spec['KE_phi'][-1:] = [En_phi]
                    spec['KE_tau'].resize((spec['KE_tau'].shape[0] + 1), axis = 0)
                    spec['KE_tau'][-1:] = [En_tau]
                    spec['KE_w'].resize((spec['KE_w'].shape[0] + 1), axis = 0)
                    spec['KE_w'][-1:] = [En_w]
                    #flux['time'].resize((flux['time'].shape[0] + 1), axis = 0)
                    #flux['time'][-1:] = [solver.sim_time]
                    #flux['KE'].resize((flux['KE'].shape[0] + 1), axis = 0)
                    #flux['KE'][-1:] = [Pi]

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
