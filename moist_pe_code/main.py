"""
Dedalus script for 3D moist primitive equations.
This script uses a Fourier basis in the horizontal directions with periodic boundary
conditions and Chebychev Polynomials in the vertical. 
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
z_basis = de.Chebyshev('z',Nz, interval=(0, Lz), dealias=3/2)
#dmain = de.Domain([x_basis, y_basis], grid_dtype=np.float64,mesh=(nnodes,16))
# Domain with nnodes commented out

domain = de.Domain([x_basis, y_basis, z_basis], grid_dtype=np.float64)

# For general use
x = domain.grid(0)
y = domain.grid(1)
z = domain.grid(2)
X,Y = np.meshgrid(x,y)
kx = domain.elements(0)
ky = domain.elements(1)
kz = domain.elements(2)
k2 = kx**2 + ky**2 + kz**2

###########
# FORCING #
###########
# Initialize fields
Fp = domain.new_field() # Phi frocing
Ft = domain.new_field() # Tau forcing

# Forcing range
#cond = (k2<=kf_max**2)&(k2>=kf_min**2)

#local_coeff_shape=Fp['c'].shape
#phase = np.random.uniform(low=-np.pi,high=np.pi,size=local_coeff_shape)
#Fp['c'][cond] = (np.cos(phase[cond])+1j*np.sin(phase[cond]))/np.sqrt(k2[cond])  # Dividing by |k| because the curl will give velocity

#local_coeff_shape=Ft['c'].shape
#phase = np.random.uniform(low=-np.pi,high=np.pi,size=local_coeff_shape)
#Ft['c'][cond] = (np.cos(phase[cond])+1j*np.sin(phase[cond]))/np.sqrt(k2[cond])

# Normalize:
#KE_op = (Fp*Fp)/(Lx*Ly) # Calculating the volume average of |Fp|^2
#KE_int_op = de.operators.integrate(KE_op,'x','y')
#E = KE_int_op.evaluate()['g'][0,0]
#Fp['g'] = Fp['g']*fp0/np.sqrt(E)

#KE_op = (Ft*Ft)/(Lx*Ly) # Calculating the volume average of |Fp|^2
#KE_int_op = de.operators.integrate(KE_op,'x','y')
#E = KE_int_op.evaluate()['g'][0,0]
#Ft['g'] = Ft['g']*ft0/np.sqrt(E)


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
problem = de.IVP(domain, variables=['u','v','w','theta','phi'])
# reduction factor
problem.parameters['r'] = r
problem.parameters['rdry'] = 1
# boundary tilt
problem.parameters['h'] = h
# beta value
problem.parameters['beta'] = beta
# drag value
problem.parameters['drag'] = drag
# epsilon
problem.parameters['eps'] = eps
# relaxation
problem.parameters['relax'] = relax
# Dissipation
problem.parameters['nu'] = nu
problem.parameters['nu_theta'] = nu_theta
problem.substitutions['nn'] = "%f" % nn # Hyperviscosity exponent
# Forcing
problem.parameters['Fp'] = Fp
problem.parameters['Ft'] = Ft
# Inverse Diffusion
problem.parameters['hnu'] = hnu
problem.substitutions['hnn'] = "%f" % hnn
# Other fields/functions
problem.parameters['h'] = h
# pi problem parameter
problem.parameters['pi'] = np.pi
# Domain size
problem.parameters['Lx'] = Lx;
problem.parameters['Ly'] = Ly;

# Operator substitutions
problem.substitutions['J(a,b)'] = "dx(a)*dy(b)-dy(a)*dx(b)"
problem.substitutions['Lap(a)'] = "dx(dx(a)) + dy(dy(a))"
problem.substitutions['Lap3D(a)'] = "dx(dx(a)) + dy(dy(a)) + dz(dz(a))"
problem.substitutions['Diff(a)'] = "Lap(Lap(a))"
problem.substitutions['Adv(a)'] = "dx(dx(dx(a))+dy(dy(a)))"
problem.substitutions['FL'] = "FractionalLaplacian"
problem.substitutions['f'] = "1"
problem.substitutions['zeta'] = "dx(v)-dy(u)"
problem.substitutions['Q'] = "(1+eps*zeta)*(1+eps*dz(theta)) - eps*eps*dz(v)*dx(theta) + eps*eps*dz(u)*dy(theta)"
problem.substitutions['theta_dot'] = "(1-thresh_r(w,r))*w*(1+eps*dz(theta))"
problem.substitutions['theta_dot_z'] = "dz(theta_dot)"
problem.substitutions['Q_dot_z'] = "eps*(1+eps*zeta)*theta_dot_z"
problem.substitutions['theta_mean'] = "-Ly/(2*pi)*cos(2*pi*y/Ly)"

# Equations of motion
problem.add_equation("eps*dt(u) + nu*Diff(u) + drag*u - f*v + dx(phi) = - eps*u*dx(u) - eps*v*dy(u) - eps*eps*w*dz(u)")
problem.add_equation("eps*dt(v) + nu*Diff(v) + drag*v + f*u + dy(phi)  = z  - eps*u*dx(v) - eps*v*dy(v) - eps*eps*w*dz(v)")
problem.add_equation("dt(theta)  + nu_theta*Diff(theta) + relax*theta - v + w = - u*dx(theta) - v*dy(theta) - eps*w*dz(theta) + (1-thresh_r(w,r))*w + eps*(1-thresh_r(w,r))*w*dz(theta)") 
problem.add_equation("dz(phi)-theta = 0")
problem.add_equation("dx(u)+dy(v)+eps*dz(w) = 0") #, condition="(nx != 0) or (ny != 0)")

# Add boundary conditions in z
problem.add_bc("left(w) = 0") #condition="(nx != 0) or (ny != 0)")
problem.add_bc("right(w) = 0",condition="(nx != 0) or (ny != 0)") #, condition="(nx != 0) or (ny != 0)")
problem.add_bc("right(phi) = 0",condition="(nx == 0) and (ny == 0)") #, condition="(nx != 0) or (ny != 0)")

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
    u = solver.state['u']
    v = solver.state['v']
    w = solver.state['w']
    phi = solver.state['phi']
    theta = solver.state['theta']

    local_coeff_shape=u['g'].shape
    u['g'] = NA*np.random.uniform(low=-1,high=1,size=local_coeff_shape)
    v['g'] = NA*np.random.uniform(low=-1,high=1,size=local_coeff_shape)
    w['g'] = NA*np.random.uniform(low=-1,high=1,size=local_coeff_shape)
    phi['g'] = NA*np.random.uniform(low=-1,high=1,size=local_coeff_shape)
    theta['g'] = NA*np.random.uniform(low=-1,high=1,size=local_coeff_shape)

    # Take away small-scale noise
    u['c'][k2>3**2] = 0.0
    v['c'][k2>3**2] = 0.0
    w['c'][k2>3**2] = 0.0
    phi['c'][k2>3**2] = 0.0
    theta['c'][k2>3**2] = 0.0

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
CFL.add_velocities(('sqrt(r)/eps','sqrt(r)/eps','sqrt(r)/eps'))

# Momentum diffusion
kcut = Nx/2. # not Nx/3 because of the 3/2 dealiasing rule instead of the 2/3 rule.
#CFL.add_frequency(nu*kcut**(2*nn))
#CFL.add_frequency(hnu*kcut**(2*hnn))

## ADD other nonlinear terms

#logger.info('1/dt restriction from visc = %f' % (nu*kcut**(2*nn)))

############
# Analysis #
############
# SNAPSHOTS
snapshots = solver.evaluator.add_file_handler('snapshots', sim_dt=tsnap_sim, max_writes=1000, mode=fh_mode)
snapshots.add_system(solver.state)
snapshots.add_task("zeta" ,name='zeta')
snapshots.add_task("Q" ,name='Q')
snapshots.add_task("u*dx(Q)" ,name='Q_adv_x')
snapshots.add_task("v*dy(Q)" ,name='Q_adv_y')
snapshots.add_task("w*dz(Q)" ,name='Q_adv_z')
snapshots.add_task("Q_dot_z" ,name='Q_dot_z')
snapshots.add_task("theta_dot" ,name='theta_dot')
snapshots.add_task("theta_dot_z" ,name='theta_dot_z')
snapshots.add_task("dz(theta)" ,name='theta_z')

## TIME SERIES
t_series = solver.evaluator.add_file_handler('time_series', iter=tseries,mode=fh_mode)

#Energy
t_series.add_task("integ(u**2 + v**2 + w**2)", name='en')

#t-equation
#t_series.add_task("integ(J(phi,tau)**2)", name='t_adv')
#t_series.add_task("integ(w**2)", name='t_w')
#t_series.add_task("integ(((1-thresh_r(w,r))*w)**2)", name='t_diab')
#t_series.add_task("integ(dx(phi)**2)", name='t_mean')
#t_series.add_task("integ((nu*FL(tau,nn))**2)", name='t_diss')

#phi-equation
#t_series.add_task("integ((J(phi,Lap(phi)))**2)", name='phi_phi')
#t_series.add_task("integ((J(tau,Lap(tau)))**2)", name='phi_tau')
#t_series.add_task("integ((h*dx(tau))**2)", name='phi_tilt')
#t_series.add_task("integ((beta*dx(phi))**2)", name='phi_beta')
#t_series.add_task("integ((Adv(tau))**2)", name='phi_adv')
#t_series.add_task("integ((nu*FL(Lap(phi),nn))**2)", name='phi_diffusion')

#tau-equation
#t_series.add_task("integ((J(tau,Lap(phi)))**2)", name='tau_tau')
#t_series.add_task("integ((J(phi,Lap(tau)))**2)", name='tau_phi')
#t_series.add_task("integ((h*dx(phi))**2)", name='tau_tilt')
#t_series.add_task("integ((beta*dx(tau))**2)", name='tau_beta')
#t_series.add_task("integ((Adv(phi))**2)", name='tau_adv')
#t_series.add_task("integ((nu*FL(Lap(tau),nn))**2)", name='tau_diffusion')
#t_series.add_task("integ((w)**2)", name='tau_w')

#dissipation
#t_series.add_task("integ(phi*nu*FL(Lap(phi),nn))",name = 'phi_diss')
#t_series.add_task("integ(phi*drag*(Lap(phi)-Lap(tau)))",name = 'phi_drag')
#t_series.add_task("integ(tau*nu*FL(Lap(tau),nn))",name = 'tau_diss')
#t_series.add_task("integ(tau*drag*(Lap(phi)-Lap(tau)))",name = 'tau_drag')

# Flow properties
flow = flow_tools.GlobalFlowProperty(solver, cadence=10)
flow.add_property("(u**2 + v**2 + w**2)/2", name='KE')
flow.add_property("(w**2)", name='KE_w')

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
            logger.info('KE = %f' %flow.volume_average('KE'))
            logger.info('KE_w = %f' %flow.volume_average('KE_w'))
        # Spectra and fluxes
        #if (solver.iteration-1) % tspec == 0:
            #logger.info('Calculating fluxes and spectra,Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
            # Load fields
            #phi  = solver.state['phi']
            #w = solver.state['w']
            #phix = x_basis.Differentiate(phi)
            #phiy = y_basis.Differentiate(phi)

            # Calculate KE spectrum
            #[modn,En_phi] = spec1d.spec1d(domain,[phix,phiy])
            #[modn,En_tau] = spec1d.spec1d(domain,[taux,tauy])
            #[modn,En_w] = spec1d.spec1d(domain,[w])
            # Calculate nonlinear term:
            #a,b,c = nonlin_term.nonlin_term(domain,[u,v,w,p])

            # Calculate KE flux:
            #[modn_f,Pi] = flux1d.flux1d_classic(domain,[u,v,w],[a,b,c],flux=True)

            # Save data:
            #if rank==0:
                #if ((solver.iteration-1)==0)&(fh_mode=='overwrite'):
                    #spec.create_dataset('time', data = [solver.sim_time],maxshape=(None,),chunks=True)
                    #spec.create_dataset('modn',data=modn)
                    #spec.create_dataset('KE_phi',data=[En_phi],maxshape=(None,np.shape(En_phi)[0]),chunks=True)
                    #spec.create_dataset('KE_tau',data=[En_tau],maxshape=(None,np.shape(En_tau)[0]),chunks=True)
                    #spec.create_dataset('KE_w',data=[En_w],maxshape=(None,np.shape(En_w)[0]),chunks=True)
                    #flux.create_dataset('time', data = [solver.sim_time],maxshape=(None,),chunks=True)
                    #flux.create_dataset('modn',data=modn_f)
                    #flux.create_dataset('KE',data=[Pi],maxshape=(None,np.shape(Pi)[0]),chunks=True)
                #else:
                    #spec['time'].resize((spec['time'].shape[0] + 1), axis = 0)
                    #spec['time'][-1:] = [solver.sim_time]
                    #spec['KE_phi'].resize((spec['KE_phi'].shape[0] + 1), axis = 0)
                    #spec['KE_phi'][-1:] = [En_phi]
                    #spec['KE_tau'].resize((spec['KE_tau'].shape[0] + 1), axis = 0)
                    #spec['KE_tau'][-1:] = [En_tau]
                    #spec['KE_w'].resize((spec['KE_w'].shape[0] + 1), axis = 0)
                    #spec['KE_w'][-1:] = [En_w]
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
