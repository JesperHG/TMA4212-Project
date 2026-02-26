"""
PDE solvers for the pure heat equation and the coupled
heat + vulcanization system.
"""

import numpy as np
from scipy.sparse.linalg import spsolve
import scipy.sparse as sparse
from src.discretization import laplacian, alpha_dot, Linear_alpha_dot
from src.params import c,gamma,T0,Tm,Ts


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

    A = laplacian(M, dx, dt)

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

    A = laplacian(M, dx, dt)

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


def solve_LinearVulcanization(M=51, N=101, T_end=15, Lambda = 1):
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

    U     = np.zeros((M, M, N))
    alpha = np.zeros((M, M, N))

    U[:, :, 0]     = -1.0
    alpha[:, :, 0]  = 0.1 

    #--------Constructing linear system----------# 
    d_star = dx**2 / dt + 4 * (1 + gamma * dx + + Lambda*(T0-Tm))
    d      = dx**2 / dt + 4 * (1 + 0.5 * gamma * dx + Lambda*(T0-Tm))

    D = sparse.diags([-1, d, -1], [-1, 0, 1], shape=(M, M), format="lil")
    D[0, 0]   = D[-1, -1] = d_star
    D[0, 1]   = D[-1, -2] = -2

    D_inner = sparse.diags([-1, 4 + dx**2 / dt + Lambda(T0-Tm), -1], [-1, 0, 1],
                           shape=(M, M), format="lil")
    D_inner[0, 0]   = D_inner[-1, -1] = d
    D_inner[0, 1]   = D_inner[-1, -2] = -2

    Amat = sparse.kron(sparse.eye(M), D_inner)
    C    = sparse.diags([-1, -1], [-M, M], shape=(M**2, M**2), format="lil")
    Amat = (Amat + C).tocsr()

    # Special cases for boundary rows
    Amat[0:M, 0:M]     = D
    Amat[0:M, M:2*M]   = 2 * Amat[0:M, M:2*M]
    Amat[-M:, -M:]      = D
    Amat[-M:, -2*M:-M]  = Amat[0:M, M:2*M]

    for i in range(N - 1):
        # Step 1: pure heat (backward Euler, 5-point stencil)
        b = (dx**2 / dt) * U[:, :, i].flatten() + Lambda*(Tm-Ts) - gamma*alpha[:,:i]
        U_star = spsolve(A, b).reshape(M, M)

        # Step 2: reaction ODE using U_star
        adot = Linear_alpha_dot(U_star, alpha[:, :, i])
        alpha[:, :, i + 1] = alpha[:, :, i] + dt * adot

        # Step 3: heat + reaction source
        b += c * dx**2 * adot.flatten()
        U[:, :, i + 1] = spsolve(A, b).reshape(M, M)

    return x, y, t, U, alpha