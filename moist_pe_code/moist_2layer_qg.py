#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""beta plane barotropic vorticity model.

Solve the barotropic vorticity equation in two dimensions

    D/Dt[ω] = 0                                                             (1)

where ω = ζ + f is absolute vorticity.  ζ is local vorticity ∇ × u and
f is global rotation.

Assuming an incompressible two-dimensional flow u = (u, v),
the streamfunction ψ = ∇ × (ψ êz) can be used to give (u,v)

    u = -∂/∂y[ψ]         v = ∂/∂x[ψ]                                        (2)

and therefore local vorticity is given by the Poisson equation

    ζ = ∆ψ                                                                  (3)

Since ∂/∂t[f] = 0 equation (1) can be written in terms of the local vorticity

        D/Dt[ζ] + u·∇f = 0
    =>  D/Dt[ζ] = -vβ                                                       (4)

using the beta-plane approximation f = f0 + βy.  This can be written entirely
in terms of the streamfunction and this is the form that will be solved
numerically.

    D/Dt[∆ψ] = -β ∂/∂x[ψ]                                                   (5)

"""

import os
import sys
import time

import numpy as np
from matplotlib import pyplot as plt
import dedalus.public as de
from dedalus.extras import flow_tools
from dedalus.tools import post

#root = logging.root
#for h in root.handlers:
#    h.setLevel("INFO")

import logging
logger = logging.getLogger(__name__)

N = 192
Lx, Ly = (1., 1.)
nx, ny = (N, N)
beta = 8.0
U = 0.0

# setup the domain
x_basis = de.Fourier('x', nx, interval=(0, Lx), dealias=3/2)
y_basis = de.Fourier('y', ny, interval=(0, Ly), dealias=3/2)
domain  = de.Domain([x_basis, y_basis], grid_dtype=np.float64)

# Stochastic forcing
def F1(*args):
    #x   = args[0].data # this is an array; we use .data to get its values
    #y   = args[1].data
    zeta   = args[0].data
    #print('size of zeta is = ', zeta.shape)
    #print('size of stochastic forcing is = ', np.random.randn(zeta.shape[0],zeta.shape[1]).shape)
    #return np.random.randn((N,N))
    #return np.random.randn()
    return 0.1*np.random.randn(zeta.shape[0],zeta.shape[1])

def SF(*args, domain=domain, F=F1):
    """This function takes arguments *args, a function F, and a domain and
    returns a Dedalus GeneralFunction that can be applied.

    """
    return de.operators.GeneralFunction(domain, layout='g', func=F1, args=args)

# now we make it parseable, so the symbol BF can be used in equations
# and the parser will know what it is.
de.operators.parseables['SF'] = SF

# Declare variables
#problem = de.IVP(domain, variables=['psi'])

problem = de.IVP(domain, variables=['phi','tau','phix','phiy','taux','tauy','w'])

# Everytime you ask for one of the expression on the left, you will get the expression on the right.
problem.substitutions['Dtau']  = " d(tau,x=2) + d(tau,y=2)"
problem.substitutions['Dphi']  = " d(phi,x=2) + d(phi,y=2)"
problem.substitutions['Dw'] = "d(w,x=2) + d(w,y=2)"

# This pattern matches for the 'thing' arguements. They don't have to be called 'thing'.
problem.substitutions['L(thing_1)']         = "  d(thing_1,x=2) + d(thing_1,y=2) "
problem.substitutions['J(thing_1,thing_2)'] = "  dx(thing_1)*dy(thing_2) - dy(thing_1)*dx(thing_2) "

# You can combine things if you want
problem.substitutions['HD(thing_1)']         = "  -D*L(L(thing_1)) "

problem.add_equation("dt(Dphi) = -J(tau,Dtau)-J(phi,Dphi)+SF(Dphi)",condition="(nx !=0) or (ny != 0)")
problem.add_equation("dt(Dtau) + w  = -J(tau,Dphi)-J(phi,Dtau)+SF(Dtau)",condition="(nx !=0) or (ny != 0)")
problem.add_equation("Dw - w = 2*J(tau,Dphi)- J(phix,taux)-J(phiy,tauy)",condition="(nx !=0) or (ny !=0)")
problem.add_equation("phix - dx(phi) = 0")
problem.add_equation("phiy - dy(phi) = 0")
problem.add_equation("taux - dx(tau) = 0")
problem.add_equation("tauy - dy(tau) = 0")
problem.add_equation("phi = 0"                         , condition="(nx == 0) and (ny == 0)")
problem.add_equation("tau = 0"                         , condition="(nx == 0) and (ny == 0)")
problem.add_equation("w = 0"                         , condition="(nx == 0) and (ny == 0)")

# Build solveru
ts = de.timesteppers.RK443
#ts = de.timesteppers.CNAB2
solver =  problem.build_solver(ts)
logger.info('Solver built')

# build the domain
x = domain.grid(0)
y = domain.grid(1)
phi = solver.state['phi']
tau = solver.state['tau']
w = solver.state['w']
phix = solver.state['phix']
phiy = solver.state['phiy']
taux = solver.state['taux']
tauy = solver.state['tauy']
tau['g'] = -y
phi['g'] = 0
w['g'] = 0;
tau.differentiate('x',out=taux)
tau.differentiate('y',out=tauy)
phi.differentiate('x',out=phix)
phi.differentiate('y',out=phiy)



# Stoptimes
solver.stop_sim_time  = np.inf
solver.stop_wall_time = np.inf
solver.stop_iteration = 300
initial_dt = 1e-3 #Lx/nx
# You can set parameters to limit the size of the timestep.
#cfl = flow_tools.CFL(solver, initial_dt=dt, cadence=10, safety=2,
#                     max_change=1.5, min_change=0.5, max_dt=10*dt)
cfl = flow_tools.CFL(solver,initial_dt,safety=0.75)
cfl.add_velocities(('phix+taux','phiy+tauy'))
dt = cfl.compute_dt()

data_dir = 'test'
#data_dir = 'active_tracer'
if domain.distributor.rank == 0:
    if not os.path.exists('{:s}/'.format(data_dir)):
         os.mkdir('{:s}/'.format(data_dir))

snapshots = solver.evaluator.add_file_handler(os.path.join(data_dir,'snapshots'), iter=10, max_writes=np.inf) #Wsim_dt=21600.,
snapshots.add_system(solver.state)
snapshots.add_task("  phi    "     ,name='phi')
#snapshots.add_task("  Dtau    "     ,name='zeta_bc')
#snapshots.add_task("  w       "     ,name='w')


analysis_tasks = [snapshots]

 # Main loop
try:
    logger.info('Starting loop')
    start_time = time.time()
    while solver.ok:
        cfldt = cfl.compute_dt()
        dt = 1e-3 # for fixed dt
        solver.step(cfldt)
        #solver.step(dt)
        if (solver.iteration-1) % 20 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e, cfldt: %e' %(solver.iteration, solver.sim_time, dt, cfldt))
        if cfldt < 1e-20:
            logger.error("cfldt very small (%e). probably stability issue" % cfldt)
            successStatus = False
            break
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_time-start_time))
    logger.info('Run time: %f cpu-hr' %((end_time-start_time)/60/60*domain.dist.comm_cart.size))


# write last checkpoint RESTART: https://groups.google.com/forum/#!searchin/dedalus-users/add_system(solver.state)%7Csort:date/dedalus-users/fUj_q91F9bo/  dcxf9NjKAgAJ
end_world_time = solver.get_world_time()
end_wall_time = end_world_time - solver.start_time
solver.evaluator.evaluate_handlers([snapshots,], timestep = dt, sim_time = solver.sim_time, world_time=end_world_time, wall_time=end_wall_time,            iteration=solver.iteration)

for task in analysis_tasks:
    logger.info(task.base_path)
    post.merge_analysis(task.base_path)
