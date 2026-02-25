import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def plotU(x, y, t, U):
    """
    Animate U(x, y, t) as a color-coded 2D plot over time.

    Parameters
    ----------
    x : 1D array of x grid points (length M+1)
    y : 1D array of y grid points (length M+1)
    t : 1D array of time points (length N+1 or N)
    U : 3D array of shape (len(x), len(y), len(t)-1) or (len(x), len(y), len(t))

    Returns
    -------
    fig, ax, anim
    """
    X, Y = np.meshgrid(x, y, indexing='ij')

    vmin = U.min()
    vmax = U.max()

    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(X, Y, U[:, :, 0], shading='auto', cmap='hot_r',
                         vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(mesh, ax=ax, label='U')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    title = ax.set_title(f't = {t[0]:.4f}')

    n_frames = U.shape[2]

    def update(frame):
        mesh.set_array(U[:, :, frame].ravel())
        title.set_text(f't = {t[frame]:.4f}')
        return mesh, title

    anim = FuncAnimation(fig, update, frames=n_frames, interval=500, blit=False)

    return fig, ax, anim


def plotAlpha(x, y, t, alpha):
    """
    Animate the degree of vulcanization alpha(x, y, t) as a color-coded 2D plot over time.

    Parameters
    ----------
    x     : 1D array of x grid points (length M+1)
    y     : 1D array of y grid points (length M+1)
    t     : 1D array of time points
    alpha : 3D array of shape (len(x), len(y), len(t)-1) with values in [0, 1]

    Returns
    -------
    fig, ax, anim
    """
    X, Y = np.meshgrid(x, y, indexing='ij')

    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(X, Y, alpha[:, :, 0], shading='auto', cmap='viridis',
                         vmin=0, vmax=1)
    cbar = fig.colorbar(mesh, ax=ax, label='Degree of vulcanization α')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    title = ax.set_title(f't = {t[0]:.4f}')

    n_frames = alpha.shape[2]

    def update(frame):
        mesh.set_array(alpha[:, :, frame].ravel())
        title.set_text(f't = {t[frame]:.4f}')
        return mesh, title

    anim = FuncAnimation(fig, update, frames=n_frames, interval=25, blit=False)

    return fig, ax, anim
