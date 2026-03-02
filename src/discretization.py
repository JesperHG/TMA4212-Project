"""
Spatial and reaction-rate discretisation helpers.
"""

import numpy as np
import scipy.sparse as sparse

from src.params import gamma, Ts, Tm, T0, E, A, R, m, n, tdiff


def sys_matrix(M, dx, dt):
    """
    Construct the discrete system matrix for the 5-point stencil with
    Robin boundary conditions.
    """
    d_star = dx**2 / dt + 4 * (1 + gamma * dx)
    d      = dx**2 / dt + 4 * (1 + 0.5 * gamma * dx)

    D = sparse.diags([-1, d, -1], [-1, 0, 1], shape=(M, M), format="lil")
    D[0, 0]   = D[-1, -1] = d_star
    D[0, 1]   = D[-1, -2] = -2

    D_inner = sparse.diags([-1, 4 + dx**2 / dt, -1], [-1, 0, 1],
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

    return Amat

def sys_matrix_3D(M, dx, dt):
    """
    Construct the discrete system matrix for the 7-point stencil in 3D with
    Robin boundary conditions.
    """
    # Inner matrix
    # Core is inner points, boundaries are side points
    inner = sparse.diags([-1, 6 + dx**2 / dt, -1], [-1, 0, 1],
                            shape=(M, M), format="lil")
    inner[0, 0] = inner[-1, -1] = 6 + dx**2 / dt + 2*dx*gamma
    inner[0, 1] = inner[-1,-2] = -2
    
    # Side matrix
    # Core is side points, boundaries are edge points
    side = sparse.diags([-1, 6 + dx**2 / dt + 2*dx*gamma, -1], [-1, 0, 1],
                            shape=(M, M), format="lil")
    side[0, 0] = side[-1,-1] = 6 + dx**2 / dt + 4*dx*gamma
    side[0, 1] = side[-1,-2] = -2

    # Edge matrix
    # Core is edge points, boundaries are corner points
    edge = sparse.diags([-1, 6 + dx**2 / dt + 4*dx*gamma, -1], [-1, 0, 1],
                            shape=(M, M), format="lil")
    edge[0, 0] = edge[-1,-1] = 6 + dx**2 / dt + 6*dx*gamma
    edge[0, 1] = edge[-1,-2] = -2
    
    # --- The hamburger model ---
    # Meat matrix
    # Core is inner matrices, boundaries are side matrices
    meat = sparse.kron(sparse.eye(M), inner) + sparse.diags(
        [-1, -1], [-M, M], shape=(M**2, M**2), format="lil")
    meat[:M,:M] = meat[-M:,-M:] = side
    meat[:M, M:2*M] = meat[-M:, -2*M:-M] = 2*meat[:M, M:2*M]

    # Bread matrix
    # Core is side matrices, boundaries are edge matrices
    bread = sparse.kron(sparse.eye(M), side) + sparse.diags(
        [-1, -1], [-M, M], shape=(M**2, M**2), format="lil")
    bread[:M,:M] = bread[-M:,-M:] = edge
    bread[:M, M:2*M] = bread[-M:, -2*M:-M] = 2*bread[:M, M:2*M]

    # Hamburger matrix (The system matrix A)
    # Core is meat matrices, boundaries are bread matrices
    hamburger = sparse.kron(sparse.eye(M), meat) + sparse.diags(
        [-1, -1], [-M**2, M**2], shape=(M**3, M**3), format="lil")
    hamburger[:M**2,:M**2] = hamburger[-M**2:,-M**2:] = bread
    hamburger[:M**2, M**2:2*M**2] = hamburger[-M**2:, -2*M**2:-M**2] = 2*hamburger[:M**2, M**2:2*M**2]

    print("----- Hamburger made -----")

    return hamburger.tocsr()


def alpha_dot(u, alpha):
    """
    Compute the rate of change of the degree of cure dα/dt (dimensionless).
    """
    limit = (Ts - Tm) / (Tm - T0)
    da   = np.zeros_like(alpha)
    mask = u > limit
    da[mask] = (A * np.exp(-E / (R * ((Tm - T0) * u[mask] + Tm)))
                * alpha[mask]**m * (1 - alpha[mask])**n)
    return da * tdiff 


def Linear_alpha_dot(u, alpha,Lambda):
    """
    Compute the rate of change of the degree of cure dα/dt (dimensionless).
    """
    limit = (Ts - Tm) / (Tm - T0)
    da   = np.zeros_like(alpha)
    mask = u > limit
    da[mask] = Lambda*(u*(Tm-T0)+Tm-Ts) - gamma * alpha
    return da * tdiff 

