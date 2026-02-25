import numpy as np
from scipy.sparse.linalg import spsolve
from System import laplacian, alphaDot
from MODEL_PARAMS import c

def solveVulcanization(M=50+1, N=100+1, T_end=15):
    '''
    Solve Heat problem with a reaction source term with the 5-point formula
    M: Number of points in each spacial direction
    N: Number of points in time
    T_end: Final dimensionless time (in units of t_diff)
    '''
    # Make the grid
    x = np.linspace(0, 1, M)
    y = np.linspace(0, 1, M) 
    t = np.linspace(0, T_end, N)
    dx = 1/(M-1)
    dt = T_end/(N-1)

    # Construct the linear system A
    A = laplacian(M, dx, dt)

    # Set up solution and source term
    U = np.zeros((M, M, N))
    alpha = np.zeros((M, M, N))

    U[:,:,0] = np.ones((M, M)) * -1
    alpha[:,:,0] = np.ones((M, M)) * 0.1

    # Solve the linear system (sparse solver) 
    # Splitting scheme
    for i in range(N-1):
        # Step 1: One step for pure heat equation (no reaction term)
        # Forward in time, 5-point in space
        b = dx**2/dt * U[:,:,i].flatten()
        U_star = spsolve(A, b).reshape(M, M) 

        # Step 2: One step for reaction term using updated solution U_star
        # Forward in time
        alphadot = alphaDot(U_star, alpha[:,:,i])
        alpha[:,:,i+1] = alpha[:,:,i] + dt*alphadot

        # Step 3: Step 1 while including the reaction term
        b += c*dx**2*alphadot.flatten()
        U[:,:,i+1] = spsolve(A, b).reshape(M, M) 
        
    return x, y, t, U, alpha