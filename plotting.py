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
    mesh = ax.pcolormesh(X, Y, U[:, :, 0], shading='auto', cmap='hot',
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

    anim = FuncAnimation(fig, update, frames=n_frames, interval=50, blit=False)

    return fig, ax, anim
