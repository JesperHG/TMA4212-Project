import numpy as np
import scipy.sparse as sparse             # Sparse matrices
from scipy.sparse.linalg import spsolve   # Linear solver for sparse matrices
from MODEL_PARAMS import c, T0, Tm, T_s, E, R, A, m, n

def laplacian(M, h, k):
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

class Vulcanization(object):
    '''
    Copy of Poisson problem class from lecture. See Poisson2D.ipynb.

    Class for the Heat problem u dot - Delta u = f on Omega, with u=g at the boundaries.
    '''
    def __init__(self, f, g):
        self.f = f   # Right hand side
        self.g = g   # Boundary values
    
    def sourceTerm(self, u, alpha): 
        limit = (T_s - Tm)/(Tm - T0) 
        da = np.zeros_like(alpha)
        mask = u > limit
        da[mask] = A*np.exp(-E/(R*((Tm-T0)*u[mask]+Tm)))*alpha[mask]**m*(1-alpha[mask])**n
        return da

    def solve(self, M=75, N=100):
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

        A = laplacian(M, h, k)

        # Solve the linear system (sparse solver)
        alpha = np.ones((M+1, M+1, N)) * 0.1

        U = np.zeros((M+1, M+1, N))
        U[1:-1,1:-1,0] = np.ones((Mi,Mi)) * -1
        U[0,0:Mi]    += self.g(xi,0.)       # lower boundary
        U[Mi-1,0:Mi] += self.g(xi,1.)       # upper boundary
        U[0:Mi,0]    += self.g(0.,yi)       # left  boundary
        U[0:Mi,Mi-1] += self.g(1.,yi)       # right boundary 


        for n in range(N-1):
            # Backward in time, 5-point in space 
            b = h**2/k*U[1:-1,1:-1,n].flatten() 
            U[1:-1,1:-1,n+1] = spsolve(A, b).reshape(Mi,Mi)  

            sourceterm = self.sourceTerm(U[1:-1,1:-1,n+1], alpha[1:-1,1:-1,n])   
            U[1:-1,1:-1,n+1] = k*c*sourceterm
            
            alpha[1:-1,1:-1,n+1] = k/2*(alpha[1:-1,1:-1,n] + k*sourceterm + alpha[1:-1,1:-1,n])
            
            U[1:-1,1:-1,n+1] += k*c*sourceterm
            
        return x, y, t, U,alpha
            # Step 1: Pure heat equation (no source term), implicit in space
            b = h**2/k * U[1:-1,1:-1,n].flatten()
            U_star = spsolve(A, b).reshape(Mi, Mi)

            # Step 2: Reaction step using updated temperature U_star
            sourceterm = self.sourceTerm(U_star, alpha[1:-1,1:-1,n])
            alpha[1:-1,1:-1,n+1] = np.clip(alpha[1:-1,1:-1,n] + k*sourceterm, 0, 1)

            # Correct U for heat released by reaction
            U[1:-1,1:-1,n+1] = U_star + k*c*sourceterm
            
        return x, y, t, U, alpha