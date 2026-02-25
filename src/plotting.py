"""
Plotting helpers – animated fields, snapshots, and analysis plots.

Plotting functions were made with help of AI.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from src.params import Tm, T0, tdiff


# ── Conversion helper ──────────────────────────────────────────────────
def U_to_T(U):
    """Convert dimensionless U to temperature T [K].  U = (T - Tm)/(Tm - T0)."""
    return U * (Tm - T0) + Tm


# ── 2-D animated temperature field ────────────────────────────────────
def plot_temperature(x, y, t, U, max_frames=60, t_end_sec=None):
    """
    Animate temperature T(x, y, t) as a colour-coded 2-D plot.

    Parameters
    ----------
    max_frames : cap on the number of animation frames
    t_end_sec  : if given, only animate up to this physical time [s]

    Returns
    -------
    fig, ax, anim
    """
    T     = U_to_T(U)
    t_sec = t * tdiff

    end_idx = T.shape[2]
    if t_end_sec is not None:
        end_idx = min(int(np.searchsorted(t_sec, t_end_sec, side='right')),
                      end_idx)

    skip          = max(1, end_idx // max_frames)
    frame_indices = np.arange(0, end_idx, skip)

    X, Y = np.meshgrid(x, y, indexing='ij')
    vmin  = T[:, :, :end_idx].min()
    vmax  = T[:, :, :end_idx].max()

    fig, ax = plt.subplots()
    mesh    = ax.pcolormesh(X, Y, T[:, :, frame_indices[0]], shading='auto',
                            cmap='coolwarm', vmin=vmin, vmax=vmax)
    fig.colorbar(mesh, ax=ax, label='Temperature T [K]')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    title = ax.set_title(f't = {t_sec[frame_indices[0]]:.1f} s')

    def update(i):
        idx = frame_indices[i]
        mesh.set_array(T[:, :, idx].ravel())
        title.set_text(f't = {t_sec[idx]:.1f} s')
        return mesh, title

    anim = FuncAnimation(fig, update, frames=len(frame_indices),
                         interval=100, blit=False)
    return fig, ax, anim


# ── 2-D animated alpha field ──────────────────────────────────────────
def plot_alpha(x, y, t, alpha, max_frames=60):
    """
    Animate the degree of vulcanisation α(x, y, t).

    Parameters
    ----------
    max_frames : cap on the number of animation frames

    Returns
    -------
    fig, ax, anim
    """
    t_sec         = t * tdiff
    n_total       = alpha.shape[2]
    skip          = max(1, n_total // max_frames)
    frame_indices = np.arange(0, n_total, skip)

    X, Y = np.meshgrid(x, y, indexing='ij')

    fig, ax = plt.subplots()
    mesh    = ax.pcolormesh(X, Y, alpha[:, :, frame_indices[0]], shading='auto',
                            cmap='Oranges', vmin=0, vmax=1)
    fig.colorbar(mesh, ax=ax, label='Degree of vulcanisation α')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    title = ax.set_title(f't = {t_sec[frame_indices[0]]:.1f} s')

    def update(i):
        idx = frame_indices[i]
        mesh.set_array(alpha[:, :, idx].ravel())
        title.set_text(f't = {t_sec[idx]:.1f} s')
        return mesh, title

    anim = FuncAnimation(fig, update, frames=len(frame_indices),
                         interval=100, blit=False)
    return fig, ax, anim


# ── Snapshots (temperature + curing side by side) ─────────────────────
def plot_snapshots(x, y, t, U, alpha, snap_times_sec=None):
    """
    Static snapshots of temperature and curing at selected physical times.

    Parameters
    ----------
    snap_times_sec : list of times in seconds; if *None* six evenly-spaced
                     snapshots are chosen automatically.

    Returns
    -------
    fig, axes
    """
    t_sec = t * tdiff
    T     = U_to_T(U)

    if snap_times_sec is None:
        indices = np.linspace(0, len(t) - 1, 6, dtype=int)
    else:
        indices = [np.argmin(np.abs(t_sec - ts)) for ts in snap_times_sec]

    X, Y = np.meshgrid(x, y, indexing='ij')
    n    = len(indices)

    fig, axes = plt.subplots(2, n, figsize=(4 * n, 7), constrained_layout=True)
    if n == 1:
        axes = axes[:, np.newaxis]

    Tmin, Tmax = T.min(), T.max()

    for j, idx in enumerate(indices):
        # Temperature
        ax = axes[0, j]
        m1 = ax.pcolormesh(X, Y, T[:, :, idx], shading='auto',
                           cmap='coolwarm', vmin=Tmin, vmax=Tmax)
        ax.set_title(f't = {t_sec[idx]:.1f} s')
        ax.set_xlabel('x');  ax.set_ylabel('y')
        fig.colorbar(m1, ax=ax, label='T [K]')

        # Curing
        ax = axes[1, j]
        m2 = ax.pcolormesh(X, Y, alpha[:, :, idx], shading='auto',
                           cmap='Oranges', vmin=0, vmax=1)
        ax.set_title(f't = {t_sec[idx]:.1f} s')
        ax.set_xlabel('x');  ax.set_ylabel('y')
        fig.colorbar(m2, ax=ax, label='α')

    fig.suptitle('Temperature (top) and Curing (bottom) snapshots', fontsize=14)
    return fig, axes


# ── Cure at centre vs time ────────────────────────────────────────────
def plot_cure_center(t, alpha):
    """
    Plot α at the domain centre over time, with a 90 % marker.

    Returns
    -------
    fig, ax
    """
    t_sec        = t * tdiff
    M            = alpha.shape[0]
    mid          = M // 2
    alpha_center = alpha[mid, mid, :]

    fig, ax = plt.subplots()
    ax.plot(t_sec, alpha_center, 'b-', linewidth=2, label='α at centre')
    ax.axhline(0.9, color='r', linestyle='--', linewidth=1, label='90 % cure')

    idx90 = np.where(alpha_center >= 0.9)[0]
    if len(idx90) > 0:
        t90 = t_sec[idx90[0]]
        ax.axvline(t90, color='gray', linestyle=':', linewidth=1)
        ax.annotate(f't₉₀ = {t90:.1f} s', xy=(t90, 0.9),
                    xytext=(t90 + 0.02 * t_sec[-1], 0.82),
                    arrowprops=dict(arrowstyle='->', color='gray'),
                    fontsize=10, color='gray')
        print(f"Time to 90 % cure at centre: {t90:.2f} s")
    else:
        print("90 % cure NOT reached at the centre within the simulation time.")

    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Degree of cure α')
    ax.set_title('Cure at centre of the domain')
    ax.set_ylim(0, 1.05)
    ax.legend()
    ax.grid(True, alpha=0.3)
    return fig, ax


# ── Maximum temperature inside the domain vs time ────────────────────
def plot_max_temp(t, U):
    """
    Plot max(T) over the domain as a function of time, with the mold
    temperature marked.  The y-axis lower limit is set near Tm so the
    overshoot region is clearly visible.

    Returns
    -------
    fig, ax
    """
    t_sec = t * tdiff
    T     = U_to_T(U)
    T_max = T.max(axis=(0, 1))

    fig, ax = plt.subplots()
    ax.plot(t_sec, T_max, 'r-', linewidth=2, label='Max T in domain')
    ax.axhline(Tm, color='k', linestyle='--', linewidth=1,
               label=f'Mold temp Tm = {Tm} K')

    overheat_mask = T_max > Tm
    if np.any(overheat_mask):
        peak_idx = np.argmax(T_max)
        ax.annotate(f'Peak = {T_max[peak_idx]:.1f} K\nt = {t_sec[peak_idx]:.1f} s',
                    xy=(t_sec[peak_idx], T_max[peak_idx]),
                    xytext=(t_sec[peak_idx] + 0.05 * t_sec[-1],
                            T_max[peak_idx] - 2),
                    arrowprops=dict(arrowstyle='->', color='gray'),
                    fontsize=10)
        print(f"Overheating detected!  Peak T = {T_max[peak_idx]:.2f} K "
              f"at t = {t_sec[peak_idx]:.2f} s")

        frame    = T[:, :, peak_idx]
        idx_flat = np.argmax(frame)
        ix, iy   = np.unravel_index(idx_flat, frame.shape)
        print(f"  Location of peak (grid indices): ix={ix}, iy={iy}")
    else:
        print("No overheating: max T never exceeds Tm.")

    y_lo = max(T_max.min(), Tm - 5)
    ax.set_ylim(y_lo, T_max.max() + 2)

    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Temperature [K]')
    ax.set_title('Maximum temperature inside the domain')
    ax.legend()
    ax.grid(True, alpha=0.3)
    return fig, ax
