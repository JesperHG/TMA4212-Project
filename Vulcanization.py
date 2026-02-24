import numpy as np
import scipy.sparse as sparse             # Sparse matrices
from scipy.sparse.linalg import spsolve   # Linear solver for sparse matrices
from MODEL_PARAMS import c, T0, Tm, T_s, E, R, A, m, n, k, h,L
from Heat import laplacian

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

    def solve(self, M=75, N=50):
        '''
        Solve Heat problem with the 5-point formula
        M: Number of grid intervals in each spacial direction
        N: Number of grid intervals in time
        '''
          # Make the grid
        x = np.linspace(0, 1, M+1)
        y = np.linspace(0, 1, M+1) 
        t = np.linspace(0, 2, N+1)
        dx = 1/M
        dt = 1/N
        Mo = M+1 # Number of total points in each direction

        # --- Constructing the linear system A  ---
        Gamma = k/(h*L)
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
        U = np.zeros((Mo, Mo, N))
        U[:,:,0] = np.ones((Mo,Mo)) * -1
        alpha = np.zeros((Mo, Mo, N)) * 0.1 
        alpha[:,:,0] = np.ones((Mo, Mo)) * 0.1

            #Time stepping 
        for n in range(N-1):
            # Step 1: Pure heat equation (no source term), implicit in space
            b = dx**2/dt * U[:,:,n].flatten()
            U_star =  U[:,:,n+1]= spsolve(A,b).reshape(Mo,Mo) 

            # Step 2: Reaction step using updated temperature U_star
            sourceterm = self.sourceTerm(U_star, alpha[:,:,n])
            alpha[:,:,n+1] = alpha[:,:,n] + dt*sourceterm

            # Correct U for heat released by reaction
            U[:,:,n+1] = U_star + dt*c*sourceterm
            
        return x, y, t, U, alpha