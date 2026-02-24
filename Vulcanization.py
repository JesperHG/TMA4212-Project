import numpy as np
import scipy.sparse as sparse             # Sparse matrices
from scipy.sparse.linalg import spsolve   # Linear solver for sparse matrices
from MODEL_PARAMS import c, T0, Tm, T_s, E, R, A, m, n, k, h, L, tdiff
from Heat import laplacian

class Vulcanization(object):
    '''
    Copy of Poisson problem class from lecture. See Poisson2D.ipynb.

    Class for the Heat problem u dot - Delta u = f on Omega, with u=g at the boundaries.
    '''
    def sourceTerm(self, u, alpha): 
        limit = (T_s - Tm)/(Tm - T0) 
        da = np.zeros_like(alpha)
        mask = u > limit
        da[mask] = A*np.exp(-E/(R*((Tm-T0)*u[mask]+Tm)))*alpha[mask]**m*(1-alpha[mask])**n
        return da*tdiff

    def solve(self, M=50, N=1000, T_end=15):
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
        Gamma = h*L/k
        d_star = dx**2/dt + 4*(1+Gamma*dx)
        d = dx**2/dt + 4+2*Gamma*dx  

        D = sparse.diags([-1,d,-1],[-1,0,1],shape=(Mo, Mo), format="lil") 
        D[0,0] = D[-1,-1] = d_star 
        D[0,1] = D[-1,-2] = -2 

        D_inner = sparse.diags([-1,4 +dx**2/dt,-1],[-1,0,1],shape=(Mo, Mo), format="lil") 
        D_inner[0,0] = D_inner[-1,-1] = d 
        D_inner[0,1] = D_inner[-1,-2] = -2 

        A = laplacian(Mo, D_inner)  
        
        # --- Implementing the special cases ---
        A[0:Mo,0:Mo] = D
        A[0:Mo,Mo:2*Mo] = 2*A[0:Mo,Mo:2*Mo] 

        A[-Mo:, -Mo:] = D  
        A[-Mo:,-2*Mo:-Mo ] = A[0:Mo,Mo:2*Mo] 

        # --- Solve the linear system (sparse solver) ---
        U = np.zeros((Mo, Mo, N+1))
        U[:,:,0] = np.ones((Mo,Mo)) * -1
        alpha = np.zeros((Mo, Mo, N+1))
        alpha[:,:,0] = np.ones((Mo, Mo)) * 0.1

        # Time stepping 
        for i in range(N):
            # Step 1: Pure heat equation (no source term), implicit in space
            b = dx**2/dt * U[:,:,i].flatten()
            U_star = spsolve(A,b).reshape(Mo,Mo) 

            # Step 2: Reaction step using updated temperature U_star
            sourceterm = self.sourceTerm(U_star, alpha[:,:,i])
            alpha[:,:,i+1] = np.clip(alpha[:,:,i] + dt*sourceterm, 0, 1)

            b += c*dx**2*sourceterm.flatten()

            # Correct U for heat released by reaction
            U[:,:,i+1] = spsolve(A,b).reshape(Mo,Mo) 
            
        return x, y, t, U, alpha