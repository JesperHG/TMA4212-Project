"""
Plotting helpers – animated fields, snapshots, and analysis plots.

Plotting functions were made with help of AI.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from src.params import Tm, T0, tdiff, m as m_exp, n as n_exp


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


# ── Temperature snapshots ─────────────────────────────────────────────
def plot_T_snapshots(x, y, t_list, U_list, labels, snap_times_sec=None):
    """
    Static snapshots of temperature fields from one or more solutions.

    Parameters
    ----------
    t_list : list of 1-D time arrays (dimensionless)
    U_list : list of 3-D U arrays, one per solution
    labels : list of row labels (e.g. ['Heat', 'Vulcanization'])
    snap_times_sec : list of snapshot times [s]; auto if *None*

    Returns
    -------
    fig, axes
    """
    n_rows = len(U_list)

    # Compute physical times and temperature arrays
    t_secs = [t * tdiff for t in t_list]
    Ts_arr = [U_to_T(U) for U in U_list]

    if snap_times_sec is None:
        snap_times_sec = np.linspace(0, min(ts[-1] for ts in t_secs), 6).tolist()

    n_cols = len(snap_times_sec)
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(3.8 * n_cols, 3.5 * n_rows),
                             constrained_layout=True, squeeze=False)

    # Global colour range across all solutions
    Tmin = min(T.min() for T in Ts_arr)
    Tmax = max(T.max() for T in Ts_arr)

    X, Y = np.meshgrid(x, y, indexing='ij')

    for i, (ts, T) in enumerate(zip(t_secs, Ts_arr)):
        for j, snap_t in enumerate(snap_times_sec):
            idx = np.argmin(np.abs(ts - snap_t))
            ax  = axes[i, j]
            m1  = ax.pcolormesh(X, Y, T[:, :, idx], shading='auto',
                               cmap='coolwarm', vmin=Tmin, vmax=Tmax)
            ax.set_title(f't = {ts[idx]:.1f} s', fontsize=9)
            ax.set_xlabel('x');  ax.set_ylabel('y')
            if j == 0:
                ax.set_ylabel(f'{labels[i]}\ny')
            fig.colorbar(m1, ax=ax, label='T [K]')

    fig.suptitle('Temperature snapshots', fontsize=14)
    return fig, axes


# ── Alpha snapshots ───────────────────────────────────────────────────
def plot_alpha_snapshots(x, y, t, alpha, snap_times_sec=None):
    """
    Static snapshots of the curing field α.

    Parameters
    ----------
    snap_times_sec : list of snapshot times [s]; auto if *None*

    Returns
    -------
    fig, axes
    """
    t_sec = t * tdiff

    if snap_times_sec is None:
        snap_times_sec = np.linspace(0, t_sec[-1], 6).tolist()

    indices = [np.argmin(np.abs(t_sec - ts)) for ts in snap_times_sec]
    n = len(indices)

    X, Y = np.meshgrid(x, y, indexing='ij')
    fig, axes = plt.subplots(1, n, figsize=(3.8 * n, 3.5),
                             constrained_layout=True, squeeze=False)

    for j, idx in enumerate(indices):
        ax = axes[0, j]
        m2 = ax.pcolormesh(X, Y, alpha[:, :, idx], shading='auto',
                           cmap='Oranges', vmin=0, vmax=1)
        ax.set_title(f't = {t_sec[idx]:.1f} s', fontsize=9)
        ax.set_xlabel('x');  ax.set_ylabel('y')
        fig.colorbar(m2, ax=ax, label='α')

    fig.suptitle('Curing (α) snapshots', fontsize=14)
    return fig, axes


# ── Cure at centre & edge vs time ─────────────────────────────────────
def _alpha_equilibrium():
    """Solve α^m = (1-α)^n for α in (0,1)."""
    from scipy.optimize import brentq
    return brentq(lambda a: a**m_exp - (1 - a)**n_exp, 0.01, 0.99)


def plot_cure_center(t, alpha):
    """
    Plot α at the domain centre and at the edge over time, with
    markers for 90 % cure and the equilibrium α^m = (1−α)^n.

    Returns
    -------
    fig, ax
    """
    t_sec        = t * tdiff
    M            = alpha.shape[0]
    mid          = M // 2
    alpha_center = alpha[mid, mid, :]
    alpha_edge   = alpha[0, mid, :]       # boundary midpoint

    alpha_eq = _alpha_equilibrium()

    fig, ax = plt.subplots()
    ax.plot(t_sec, alpha_center, 'b-', linewidth=2, label='α at centre')
    ax.plot(t_sec, alpha_edge,   'g-', linewidth=2, label='α at edge')
    ax.axhline(0.9, color='r', linestyle='--', linewidth=1, label='90 % cure')
    ax.axhline(alpha_eq, color='purple', linestyle=':', linewidth=1,
               label=f'α^m = (1−α)^n  (α ≈ {alpha_eq:.3f})')

    # Mark 90 % time at centre
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
    ax.set_title('Cure at centre and edge of the domain')
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    return fig, ax


# ── Maximum temperature inside the domain vs time ────────────────────
def plot_max_temp(t, U, alpha=None):
    """
    Plot max(T) over the domain as a function of time, with the mold
    temperature marked.  Optionally overlay the centre cure α on a
    secondary y-axis to show correlation.

    Parameters
    ----------
    alpha : 3-D array, optional.  If given, centre cure is plotted on a
            twin y-axis.

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
    ax.grid(True, alpha=0.3)

    # Overlay centre cure on twin axis
    if alpha is not None:
        M_    = alpha.shape[0]
        mid   = M_ // 2
        ax2   = ax.twinx()
        ax2.plot(t_sec, alpha[mid, mid, :], 'b-', linewidth=1.5,
                 label='α at centre')
        ax2.set_ylabel('Degree of cure α', color='b')
        ax2.set_ylim(0, 1.05)
        ax2.tick_params(axis='y', labelcolor='b')
        # Merge legends from both axes
        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, fontsize=8,
                  loc='center right')
    else:
        ax.legend(fontsize=8)

    return fig, ax
