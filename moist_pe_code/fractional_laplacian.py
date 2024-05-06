"""Fourier fractional Laplacian operator."""

import numpy as np
import dedalus.public as de

from dedalus.core.basis import Fourier
from dedalus.core.field import Operand
from dedalus.core.operators import Separable, FutureField
from dedalus.tools.array import reshape_vector


class FractionalLaplacian(Separable, FutureField):
    """
    Fourier fractional Laplacian operator: (-Δ)**s

    Parameters
    ----------
    arg : field object
        Field argument
    s : float
        Laplacian power

    Notes
    -----
    The fractional Laplacian is defined as (-Δ)**s, with a corresponding
    Fourier symbol |k|**(2s).  The standard Laplacian is recovered, up to an
    overall negative sign, with s=1.

    See https://www.ma.utexas.edu/mediawiki/index.php/Fractional_Laplacian

    """

    def __new__(cls, arg0, *args, **kw):
        # Cast to operand
        arg0 = Operand.cast(arg0)
        # Check all bases are Fourier
        for basis in arg0.domain.bases:
            if not isinstance(basis, Fourier):
                raise NotImplementedError("Operator only implemented for full-Fourier domains. ")
        # Check for scalars
        if arg0.domain.dim == 0:
            return 0
        else:
            return object.__new__(cls)

    def __init__(self, arg, s, **kw):
        arg = Operand.cast(arg)
        super().__init__(arg, **kw)
        self.kw = {'s': s}
        self.s = s
        self.name = 'Lap[%s]' % self.s
        self.axis = None
        # Build operator symbol array
        slices = self.domain.dist.coeff_layout.slices(self.domain.dealias)
        local_wavenumbers = [self.domain.elements(axis) for axis in range(self.domain.dim)]
        local_k2 = np.sum(ki**2 for ki in local_wavenumbers)
        local_k2_mod = local_k2.copy()
        local_k2_mod[local_k2 == 0] = 1
        self.local_symbols = local_k2_mod ** s
        self.local_symbols[local_k2 == 0] = 0

    def meta_constant(self, axis):
        # Preserve constancy
        return self.args[0].meta[axis]['constant']

    def check_conditions(self):
        arg0, = self.args
        # Must be in coeff layout
        is_coeff = not np.any(arg0.layout.grid_space)
        return is_coeff

    def operator_form(self, index):
        # Get local index, special casing zero mode for NCC preconstruction pass
        if any(index):
            local_index = index - self.domain.dist.coeff_layout.start(scales=None)
        else:
            local_index = index
        return self.local_symbols[tuple(local_index)]

    def operate(self, out):
        arg0, = self.args
        # Require coeff layout
        arg0.require_coeff_space()
        out.layout = arg0.layout
        # Apply symbol array to coefficients
        np.multiply(arg0.data, self.local_symbols, out=out.data)


def test_eval(Nx, Ny, Lx, Ly, s):
    """
    Test fractional Laplacian evaluation and compare to analytical result.

    Parameters
    ----------
    Nx, Ny : int
        x and y resolutions
    Lx, Ly : float
        x and y box lengths
    mode : tuple of ints
        Mode number to seed and test
    s : float
        Fractional Laplacian power

    """

    # Bases and domain
    xbasis = de.Fourier('x', Nx, interval=[0, Lx])
    ybasis = de.Fourier('y', Ny, interval=[0, Ly])
    domain = de.Domain([xbasis, ybasis], grid_dtype=np.float64)

    # Input
    u = domain.new_field()
    u['g'] = np.random.randn(*u['g'].shape)

    # True solution
    kx = domain.elements(0)
    ky = domain.elements(1)
    k2 = kx**2 + ky**2
    k2_mod = k2.copy()
    k2_mod[k2 == 0] = 1
    w_true = domain.new_field()
    w_true['c'] = u['c'] * k2_mod**s
    w_true['c'][k2 == 0] = 0

    # Test operator
    w = FractionalLaplacian(u, s)
    print('Test eval:', np.allclose(w['c'], w_true['c']))


def test_ivp(Nx, Ny, Lx, Ly, s, dt, iters):
    """
    Test fractional Laplacian on 2D heat equation and compare to analytical result.

        dt(u) + FL(u, s) = 0

    Parameters
    ----------
    Nx, Ny : int
        x and y resolutions
    Lx, Ly : float
        x and y box lengths
    s : float
        Fractional Laplacian power
    dt : float
        Timestep
    iters : int
        Iterations

    """

    # Add operator to namespace
    de.operators.parseables['FractionalLaplacian'] = FractionalLaplacian

    # Bases and domain
    xbasis = de.Fourier('x', Nx, interval=[0, Lx])
    ybasis = de.Fourier('y', Ny, interval=[0, Ly])
    domain = de.Domain([xbasis, ybasis], grid_dtype=np.float64)

    # Problem
    problem = de.IVP(domain, ['u'])
    problem.substitutions['s'] = "%f" % s
    problem.substitutions['FL'] = "FractionalLaplacian"
    problem.add_equation("dt(u) + 2*FL(u, s) = u*u - u*u + FL(u,s)") # Add and subtract nonlinear terms to force transforms in RHS evaluation

    # Solver
    solver = problem.build_solver(de.timesteppers.RK222)
    solver.stop_iteration = iters

    # Initial conditions
    u = solver.state['u']
    u['g'] = np.random.randn(*u['g'].shape)
    uc0 = u['c'].copy()

    # True solution
    kx = domain.elements(0)
    ky = domain.elements(1)
    k2 = kx**2 + ky**2
    k2_mod = k2.copy()
    k2_mod[k2 == 0] = 1
    uc_true = u['c'] * np.exp(-k2_mod**s * dt * iters)
    uc_true[k2 == 0] = u['c'][k2 == 0]

    # Test problem
    while solver.ok:
        solver.step(dt)
    print('Test IVP:', np.allclose(u['c'], uc_true, atol=1e-6))


if __name__ == '__main__':

    Nx = 32
    Ny = 32
    Lx = 5
    Ly = 2
    s = -1
    dt = 1e-4
    iters = 100

    test_eval(Nx, Ny, Lx, Ly, s)
    test_ivp(Nx, Ny, Lx, Ly, s, dt, iters)

