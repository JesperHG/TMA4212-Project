import numpy as np
from scipy.sparse.linalg import spsolve
from System import laplacian

def solveHeat(M=50+1, N=100+1, T_end=1):
    '''
    Solve Heat problem with the 5-point formula
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

    # Set up solution
    U = np.zeros((M, M, N))
    U[:,:,0] = np.ones((M, M)) * -1

    # Solve the linear system (sparse solver)
    for i in range(N-1):
        # Backward in time, 5-point in space
        b = dx**2/dt * U[:,:,i].flatten()
        U[:,:,i+1]= spsolve(A, b).reshape(M, M) 

    return x, y, t, U