import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from MODEL_PARAMS import k, h, rho, c_p, H_r, A, E, R, m, n, T_s, L, T0, Tm
from Heat import Heat
from Vulcanization import Vulcanization

def plot_solution(x, y, t, U, text='The solution', fp=1):
    """
    This function was made with AI assistance.

    Animate 2D solution u(x,y,t) with x,y on axes, color-coded by U, animated over t.
    
    Parameters:
    x, y: 1D arrays of spatial coordinates
    t: 1D array of time values
    U: 3D array of shape (len(x), len(y), len(t))
    text: title for the plot
    fp: figure number
    """
    fig, ax = plt.subplots(figsize=(8, 6), num=fp)
    
    # Create meshgrid for x, y
    X, Y = np.meshgrid(x, y, indexing='ij') 

    u_max = np.max(U)
    u_min = np.min(U)
    
    # Create initial contourf plot
    contour = ax.contourf(X, Y, U[:, :, 0], levels=20, cmap='coolwarm', vmin=u_min, vmax=u_max)
    cbar = fig.colorbar(contour, ax=ax, label='U')
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'{text} (t = {t[0]:.4f})')
    
    def animate(frame):
        ax.clear()
        contour = ax.contourf(X, Y, U[:, :, frame], levels=20, cmap='coolwarm', vmin=u_min, vmax=u_max)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(f'{text} (t = {t[frame]:.4f})')
        return contour,
    
    anim = FuncAnimation(fig, animate, frames=len(t)-1, interval=100, blit=False, repeat=True)
    
    return fig, ax, anim



f = lambda x, y: 0

pde = Vulcanization(f) 
pde2 = Heat(f)

x, y, t, U = pde2.solve() 

# Create animated plot
fig, ax, anim = plot_solution(x, y, t, U)
plt.show() 

