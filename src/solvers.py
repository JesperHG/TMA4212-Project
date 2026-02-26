"""
PDE solvers for the pure heat equation and the coupled
heat + vulcanization system.
"""

import numpy as np
from scipy.sparse.linalg import spsolve

from src.discretization import sys_matrix, alpha_dot
from src.params import c


def solve_heat(M=51, N=101, T_end=1):
    """
    Solve the pure heat equation (no reaction term) with the 5-point formula.

    Parameters
    ----------
    M     : number of grid points in each spatial direction
    N     : number of time steps
    T_end : final dimensionless time (in units of t_diff)

    Returns
    -------
    x, y, t, U
    """
    x  = np.linspace(0, 1, M)
    y  = np.linspace(0, 1, M)
    t  = np.linspace(0, T_end, N)
    dx = 1 / (M - 1)
    dt = T_end / (N - 1)

    A = sys_matrix(M, dx, dt)

    U = np.zeros((M, M, N))
    U[:, :, 0] = -1.0

    for i in range(N - 1):
        b = (dx**2 / dt) * U[:, :, i].flatten()
        U[:, :, i + 1] = spsolve(A, b).reshape(M, M)

    return x, y, t, U


def solve_vulcanization(M=51, N=101, T_end=15):
    """
    Solve the coupled heat + vulcanization system using operator splitting.

    Parameters
    ----------
    M     : number of grid points in each spatial direction
    N     : number of time steps
    T_end : final dimensionless time (in units of t_diff)

    Returns
    -------
    x, y, t, U, alpha
    """
    x  = np.linspace(0, 1, M)
    y  = np.linspace(0, 1, M)
    t  = np.linspace(0, T_end, N)
    dx = 1 / (M - 1)
    dt = T_end / (N - 1)

    A = sys_matrix(M, dx, dt)

    U     = np.zeros((M, M, N))
    alpha = np.zeros((M, M, N))

    U[:, :, 0]     = -1.0
    alpha[:, :, 0]  = 0.1

    for i in range(N - 1):
        # Step 1: pure heat (backward Euler, 5-point stencil)
        b      = (dx**2 / dt) * U[:, :, i].flatten()
        U_star = spsolve(A, b).reshape(M, M)

        # Step 2: reaction ODE using U_star
        adot = alpha_dot(U_star, alpha[:, :, i])
        alpha[:, :, i + 1] = alpha[:, :, i] + dt * adot

        # Step 3: heat + reaction source
        b += c * dx**2 * adot.flatten()
        U[:, :, i + 1] = spsolve(A, b).reshape(M, M)

    return x, y, t, U, alpha
