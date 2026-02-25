import scipy.sparse as sparse
import numpy as np
from MODEL_PARAMS import gamma, Ts, Tm, T0, E, A, R, m, n, tdiff

def laplacian(M, dx, dt):
    '''
    Construct the Discrete Laplacian
    '''
    d_star = dx**2/dt + 4*(1 + gamma*dx)
    d = dx**2/dt + 4*(1 + 0.5*gamma*dx)

    D = sparse.diags([-1,d,-1], [-1,0,1], shape=(M,M), format="lil") 
    D[0,0] = D[-1,-1] = d_star 
    D[0,1] = D[-1,-2] = -2 

    D_inner = sparse.diags([-1,4 + dx**2/dt,-1], [-1,0,1], shape=(M, M), format="lil") 
    D_inner[0,0] = D_inner[-1,-1] = d 
    D_inner[0,1] = D_inner[-1,-2] = -2 

    Amat = sparse.kron(sparse.eye(M), D_inner)
    C = sparse.diags([-1,-1], [-M,M], shape=(M**2, M**2), format="lil")
    Amat = (Amat+C).tocsr() # Konverter til csr-format (necessary for spsolve) 

    # Implementing the special cases
    Amat[0:M, 0:M] = D
    Amat[0:M, M:2*M] = 2*Amat[0:M, M:2*M] 

    Amat[-M:, -M:] = D  
    Amat[-M:, -2*M:-M] = Amat[0:M, M:2*M] 

    return Amat

def alphaDot(u, alpha): 
    limit = (Ts - Tm)/(Tm - T0) 
    da = np.zeros_like(alpha)
    mask = u > limit
    da[mask] = A*np.exp(-E/(R*((Tm-T0)*u[mask]+Tm)))*alpha[mask]**m*(1-alpha[mask])**n
    return da*tdiff