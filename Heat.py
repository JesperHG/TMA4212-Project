import numpy as np
import scipy.sparse as sparse             # Sparse matrices
from scipy.sparse.linalg import spsolve
from MODEL_PARAMS import L,h,k   

def laplacian(Mo, D):
    '''
    Construct the Discrete Laplacian
    '''
    # Construct a sparse A-matrix
    A = sparse.kron(sparse.eye(Mo), D)
    C = sparse.diags([-1,-1],[-Mo,Mo],shape=(Mo**2, Mo**2), format="lil")
    A = (A+C).tocsr() # Konverter til csr-format (necessary for spsolve) 
    return A  

class Heat(object):
    '''
    Copy of Poisson problem class from lecture. See Poisson2D.ipynb.

    Class for the Heat problem u dot - Delta u = f on Omega, with u=g at the boundaries.
    '''

    def solve(self, M=50, N=100, T_end=1):
        '''
        Solve Heat problem with the 5-point formula
        M: Number of grid intervals in each spacial direction
        N: Number of grid intervals in time
        T_end: Final dimensionless time (in units of t_diff)
        '''

        # Make the grid
        x = np.linspace(0, 1, M+1)
        y = np.linspace(0, 1, M+1) 
        t = np.linspace(0, T_end, N+1)
        dx = 1/M
        dt = T_end/N
        Mo = M+1 # Number of total points in each direction

        # --- Constructing the linear system A  ---
        gamma = h*L/k
        d_star = dx**2/dt + 4*(1+gamma*dx)
        d = dx**2/dt + 4+2*gamma*dx  

        D = sparse.diags([-1,d,-1],[-1,0,1],shape=(Mo, Mo), format="lil") 
        D[0,0] = D[-1,-1] = d_star 
        D[0,1] = D[-1,-2] = -2 

        D_inner = sparse.diags([-1,4 +dx**2/dt,-1],[-1,0,1],shape=(Mo, Mo), format="lil") 
        D_inner[0,0] = D_inner[-1,-1] = d 
        D_inner[0,1] = D_inner[-1,-2] = -2 

        A = laplacian(Mo, D_inner)  
        
        # --- Implementing the special cases ---
        A[0:Mo, 0:Mo] = D
        A[0:Mo, Mo:2*Mo] = 2*A[0:Mo, Mo:2*Mo] 

        A[-Mo:, -Mo:] = D  
        A[-Mo:, -2*Mo:-Mo] = A[0:Mo, Mo:2*Mo] 

        # --- Solve the linear system (sparse solver) ---
        U = np.zeros((Mo, Mo, N+1))
        U[:,:,0] = np.ones((Mo, Mo)) * -1

        for i in range(N):
            # Backward in time, 5-point in space
            U[:,:,i+1]= spsolve(A, dx**2/dt*U[:,:,i].flatten()).reshape(Mo, Mo) 

        return x, y, t, U