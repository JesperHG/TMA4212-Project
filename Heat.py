import numpy as np
import scipy.sparse as sparse             # Sparse matrices
from scipy.sparse.linalg import spsolve
from MODEL_PARAMS import L   # Linear solver for sparse matrices

def laplacian(M, D):
    '''
    Construct the Discrete Laplacian
    '''
                       

    # Construct a sparse A-matrix
    A = sparse.kron(sparse.eye(M), D)
    C = sparse.diags([-1,-1],[-M,M],shape=(M**2, M**2), format="lil")
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
        M = M + 1

        # Inner points
        xi = x[1:-1]       
        yi = y[1:-1] 
        Mi = M-1          # Number of inner points in each direction

        #Constructing the linear system A  
        d_star = h**2/k +4*(1+k)
        d = h**2/k + 4 + 2*k  

        D = sparse.diags([-1,d,-1],[-1,0,1],shape=(M, M), format="lil") 
        D[0,0] = D[-1,-1] = d_star 
        D[0,1] = D[-1,-2] = -2 

        D_inner = sparse.diags([-1,4 +h**2/k,-1],[-1,0,1],shape=(M, M), format="lil") 
        D_inner[0,0] = D_inner[-1,-1] = d 
        D_inner[0,1] = D_inner[-1,-2] = -2 


        A = laplacian(M,D_inner)  
        
        #Implementing the special cases 
        A[0:M,0:M] = D
        A[0:M,M:2*M] = 2*A[0:M,M:2*M] 

        A[-M:, -M:] = D  
        A[-M:,-2*M:-M ] = A[0:M,M:2*M] 





        # Solve the linear system (sparse solver)
        U = np.zeros((M, M, N))
        U[:,:,0] = np.ones((M,M)) * -1
      

        for n in range(N-1):
            #Interior points
                # Backward in time, 5-point in space
            U[:,:,n+1]= spsolve(A, h**2/k*U[:,:,n].flatten()).reshape(M,M) 

        
            
        return x, y, t, U