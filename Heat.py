import numpy as np
import scipy.sparse as sparse             # Sparse matrices
from scipy.sparse.linalg import spsolve   # Linear solver for sparse matrices

def laplacian(M,h,k):
    '''
    Construct the Discrete Laplacian
    '''
    Mi = M-1       # Number of inner points in each direction
    Mi2 = Mi**2    # Number of inner points in total

    # Construct a sparse A-matrix
    B = sparse.diags([-1,4+h**2/k,-1],[-1,0,1],shape=(Mi, Mi), format="lil")
    A = sparse.kron(sparse.eye(Mi), B)
    C = sparse.diags([-1,-1],[-Mi,Mi],shape=(Mi2, Mi2), format="lil")
    A = (A+C).tocsr() # Konverter til csr-format (necessary for spsolve) 
    return A

class Heat(object):
    '''
    Copy of Poisson problem class from lecture. See Poisson2D.ipynb.

    Class for the Heat problem u dot - Delta u = f on Omega, with u=g at the boundaries.
    '''
    def __init__(self, f, g):
        self.f = f   # Right hand side
        self.g = g   # Boundary values

    def solve(self, M=10, N=50):
        '''
        Solve Heat problem with the 5-point formula
        M: Number of grid intervals in each spacial direction
        N: Number of grid intervals in time
        '''

        # Make the grid
        x = np.linspace(0, 1, M+1)
        y = np.linspace(0, 1, M+1) 
        t = np.linspace(0, 1, N+1)
        h = 1/M
        k = 1/N

        # Inner points
        xi = x[1:-1]       
        yi = y[1:-1] 
        Mi = M-1          # Number of inner points in each direction

        A = laplacian(M, k, h)

        # Solve the linear system (sparse solver)
        U = np.zeros((M+1, M+1, N))
        U[1:-1,1:-1,0] = np.ones((Mi,Mi)) * -1
        U[0,0:Mi]    += self.g(xi,0.)       # lower boundary
        U[Mi-1,0:Mi] += self.g(xi,1.)       # upper boundary
        U[0:Mi,0]    += self.g(0.,yi)       # left  boundary
        U[0:Mi,Mi-1] += self.g(1.,yi)       # right boundary

        for n in range(N-1):
            # Backward in time, 5-point in space
            U[1:-1,1:-1,n+1] = spsolve(A, h**2/k*U[1:-1,1:-1,n].flatten()).reshape(Mi,Mi)
            
        return x, y, t, U