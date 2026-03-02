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


# ═══════════════════════════════════════════════════════════════════════
#  3-D  PLOTTING  HELPERS
# ═══════════════════════════════════════════════════════════════════════

# ── Animated z-slice temperature ──────────────────────────────────────
def plot_temperature_3D_slices(x, y, z, t, U, z_fracs=(0.0, 0.25, 0.5),
                               max_frames=60, t_end_sec=None):
    """
    Animate T(x, y, z_i, t) side-by-side for several z-slices.

    Parameters
    ----------
    z_fracs    : fractional z positions to slice (0 = boundary, 0.5 = centre)
    max_frames : cap on animation frames
    t_end_sec  : physical-time upper limit [s]

    Returns
    -------
    fig, axes, anim
    """
    T     = U_to_T(U)
    t_sec = t * tdiff
    M     = U.shape[0]

    z_indices = [min(int(f * (M - 1)), M - 1) for f in z_fracs]

    end_idx = T.shape[3]
    if t_end_sec is not None:
        end_idx = min(int(np.searchsorted(t_sec, t_end_sec, side='right')),
                      end_idx)

    skip          = max(1, end_idx // max_frames)
    frame_indices = np.arange(0, end_idx, skip)

    X, Y = np.meshgrid(x, y, indexing='ij')

    vmin = T[:, :, z_indices, :end_idx].min()
    vmax = T[:, :, z_indices, :end_idx].max()

    n_slices = len(z_indices)
    fig, axes = plt.subplots(1, n_slices, figsize=(5 * n_slices, 4.5),
                             constrained_layout=True)
    if n_slices == 1:
        axes = [axes]

    meshes = []
    for ax, zi, zf in zip(axes, z_indices, z_fracs):
        m1 = ax.pcolormesh(X, Y, T[:, :, zi, frame_indices[0]],
                           shading='auto', cmap='coolwarm', vmin=vmin, vmax=vmax)
        fig.colorbar(m1, ax=ax, label='T [K]')
        ax.set_xlabel('x'); ax.set_ylabel('y')
        ax.set_title(f'z = {z[zi]:.2f}  (z/{z[-1]:.0f} = {zf:.2f})')
        meshes.append(m1)

    suptitle = fig.suptitle(f't = {t_sec[frame_indices[0]]:.1f} s', fontsize=13)

    def update(i):
        idx = frame_indices[i]
        for m1, zi in zip(meshes, z_indices):
            m1.set_array(T[:, :, zi, idx].ravel())
        suptitle.set_text(f't = {t_sec[idx]:.1f} s')
        return meshes + [suptitle]

    anim = FuncAnimation(fig, update, frames=len(frame_indices),
                         interval=100, blit=False)
    return fig, axes, anim


# ── Animated z-slice alpha ────────────────────────────────────────────
def plot_alpha_3D_slices(x, y, z, t, alpha, z_fracs=(0.0, 0.25, 0.5),
                         max_frames=60):
    """
    Animate α(x, y, z_i, t) side-by-side for several z-slices.

    Returns
    -------
    fig, axes, anim
    """
    t_sec     = t * tdiff
    M         = alpha.shape[0]
    z_indices = [min(int(f * (M - 1)), M - 1) for f in z_fracs]

    n_total       = alpha.shape[3]
    skip          = max(1, n_total // max_frames)
    frame_indices = np.arange(0, n_total, skip)

    X, Y = np.meshgrid(x, y, indexing='ij')
    n_slices = len(z_indices)
    fig, axes = plt.subplots(1, n_slices, figsize=(5 * n_slices, 4.5),
                             constrained_layout=True)
    if n_slices == 1:
        axes = [axes]

    meshes = []
    for ax, zi, zf in zip(axes, z_indices, z_fracs):
        m1 = ax.pcolormesh(X, Y, alpha[:, :, zi, frame_indices[0]],
                           shading='auto', cmap='Oranges', vmin=0, vmax=1)
        fig.colorbar(m1, ax=ax, label='α')
        ax.set_xlabel('x'); ax.set_ylabel('y')
        ax.set_title(f'z = {z[zi]:.2f}')
        meshes.append(m1)

    suptitle = fig.suptitle(f't = {t_sec[frame_indices[0]]:.1f} s', fontsize=13)

    def update(i):
        idx = frame_indices[i]
        for m1, zi in zip(meshes, z_indices):
            m1.set_array(alpha[:, :, zi, idx].ravel())
        suptitle.set_text(f't = {t_sec[idx]:.1f} s')
        return meshes + [suptitle]

    anim = FuncAnimation(fig, update, frames=len(frame_indices),
                         interval=100, blit=False)
    return fig, axes, anim


# ── 3-D Temperature snapshots (z-slices × time) ──────────────────────
def plot_T_snapshots_3D(x, y, z, t, U,
                        z_fracs=(0.0, 0.25, 0.5),
                        snap_times_sec=None):
    """
    Static snapshot grid: rows = z-slices, columns = time snapshots.

    Returns
    -------
    fig, axes
    """
    T     = U_to_T(U)
    t_sec = t * tdiff
    M     = U.shape[0]

    z_indices = [min(int(f * (M - 1)), M - 1) for f in z_fracs]

    if snap_times_sec is None:
        snap_times_sec = np.linspace(0, t_sec[-1], 6).tolist()

    n_rows = len(z_indices)
    n_cols = len(snap_times_sec)

    Tmin = T[:, :, z_indices, :].min()
    Tmax = T[:, :, z_indices, :].max()

    X, Y = np.meshgrid(x, y, indexing='ij')

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(3.8 * n_cols, 3.5 * n_rows),
                             constrained_layout=True, squeeze=False)

    for i, (zi, zf) in enumerate(zip(z_indices, z_fracs)):
        for j, snap_t in enumerate(snap_times_sec):
            tidx = np.argmin(np.abs(t_sec - snap_t))
            ax   = axes[i, j]
            m1   = ax.pcolormesh(X, Y, T[:, :, zi, tidx], shading='auto',
                                 cmap='coolwarm', vmin=Tmin, vmax=Tmax)
            ax.set_title(f't = {t_sec[tidx]:.1f} s', fontsize=9)
            ax.set_xlabel('x')
            if j == 0:
                ax.set_ylabel(f'z = {z[zi]:.2f}\ny')
            else:
                ax.set_ylabel('y')
            fig.colorbar(m1, ax=ax, label='T [K]')

    fig.suptitle('3-D Temperature snapshots (z-slices)', fontsize=14)
    return fig, axes


# ── 3-D Alpha snapshots (z-slices × time) ────────────────────────────
def plot_alpha_snapshots_3D(x, y, z, t, alpha,
                            z_fracs=(0.0, 0.25, 0.5),
                            snap_times_sec=None):
    """
    Static snapshot grid for α: rows = z-slices, columns = time snapshots.

    Returns
    -------
    fig, axes
    """
    t_sec = t * tdiff
    M     = alpha.shape[0]

    z_indices = [min(int(f * (M - 1)), M - 1) for f in z_fracs]

    if snap_times_sec is None:
        snap_times_sec = np.linspace(0, t_sec[-1], 6).tolist()

    n_rows = len(z_indices)
    n_cols = len(snap_times_sec)

    X, Y = np.meshgrid(x, y, indexing='ij')

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(3.8 * n_cols, 3.5 * n_rows),
                             constrained_layout=True, squeeze=False)

    for i, (zi, zf) in enumerate(zip(z_indices, z_fracs)):
        for j, snap_t in enumerate(snap_times_sec):
            tidx = np.argmin(np.abs(t_sec - snap_t))
            ax   = axes[i, j]
            m1   = ax.pcolormesh(X, Y, alpha[:, :, zi, tidx], shading='auto',
                                 cmap='Oranges', vmin=0, vmax=1)
            ax.set_title(f't = {t_sec[tidx]:.1f} s', fontsize=9)
            ax.set_xlabel('x')
            if j == 0:
                ax.set_ylabel(f'z = {z[zi]:.2f}\ny')
            else:
                ax.set_ylabel('y')
            fig.colorbar(m1, ax=ax, label='α')

    fig.suptitle('3-D Curing (α) snapshots (z-slices)', fontsize=14)
    return fig, axes


# ── Centre-line profile along z ───────────────────────────────────────
def plot_centerline_z(z, t, U, alpha=None,
                      snap_times_T_sec=None, snap_times_alpha_sec=None):
    """
    Plot temperature (and optionally α) along the centre-line
    x = y = 0.5 as a function of z for several time snapshots.

    Parameters
    ----------
    snap_times_T_sec     : snapshot times [s] for the temperature panel
    snap_times_alpha_sec : snapshot times [s] for the α panel (defaults to
                           *snap_times_T_sec* when not given)

    Returns
    -------
    fig, axes   (one or two subplots)
    """
    T     = U_to_T(U)
    t_sec = t * tdiff
    M     = U.shape[0]
    mid   = M // 2

    if snap_times_T_sec is None:
        snap_times_T_sec = np.linspace(0, t_sec[-1], 8).tolist()
    if snap_times_alpha_sec is None:
        snap_times_alpha_sec = snap_times_T_sec

    T_indices = [np.argmin(np.abs(t_sec - ts)) for ts in snap_times_T_sec]
    a_indices = [np.argmin(np.abs(t_sec - ts)) for ts in snap_times_alpha_sec]

    n_axes = 2 if alpha is not None else 1
    fig, axes = plt.subplots(1, n_axes, figsize=(6 * n_axes, 4.5),
                             constrained_layout=True, squeeze=False)

    ax_T = axes[0, 0]
    for idx in T_indices:
        ax_T.plot(z, T[mid, mid, :, idx], label=f't = {t_sec[idx]:.0f} s')
    ax_T.set_xlabel('z')
    ax_T.set_ylabel('Temperature [K]')
    ax_T.set_title('T along centre-line (x = y = 0.5)')
    ax_T.legend(fontsize=7, ncol=2)
    ax_T.grid(True, alpha=0.3)

    if alpha is not None:
        ax_a = axes[0, 1]
        for idx in a_indices:
            ax_a.plot(z, alpha[mid, mid, :, idx],
                      label=f't = {t_sec[idx]:.0f} s')
        ax_a.set_xlabel('z')
        ax_a.set_ylabel('Degree of cure α')
        ax_a.set_title('α along centre-line (x = y = 0.5)')
        ax_a.set_ylim(0, 1.05)
        ax_a.legend(fontsize=7, ncol=2)
        ax_a.grid(True, alpha=0.3)

    return fig, axes


# ── 3-D isosurface snapshots ─────────────────────────────────────────
def plot_isosurface_3D(x, y, z, t, U, alpha,
                       alpha_levels=None,
                       snap_times_sec_per_level=None):
    """
    3×3 grid of 3-D scatter 'isosurface' plots.

    Rows  = different α levels (below, at, and above the equilibrium).
    Columns = time snapshots (different per level, since the cure front
    progresses at different speeds for different thresholds).

    Parameters
    ----------
    alpha_levels : list of 3 floats  (default: eq−0.15, eq, eq+0.15)
    snap_times_sec_per_level : list of 3 lists, each with 3 times [s]

    Returns
    -------
    fig
    """
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    T     = U_to_T(U)
    t_sec = t * tdiff

    alpha_eq = _alpha_equilibrium()

    if alpha_levels is None:
        alpha_levels = [alpha_eq - 0.15, alpha_eq, alpha_eq + 0.15]
    if snap_times_sec_per_level is None:
        snap_times_sec_per_level = [
            [10000, 20000, 40000],   # below equilibrium – reached early
            [20000, 40000, 60000],   # at equilibrium
            [40000, 60000, 80000],   # above equilibrium – reached late
        ]

    n_rows = len(alpha_levels)
    n_cols = max(len(ts) for ts in snap_times_sec_per_level)

    X3, Y3, Z3 = np.meshgrid(x, y, z, indexing='ij')
    Tmin_global = Tm - 5
    Tmax_global = T.max()

    level_tags = []
    for i, alev in enumerate(alpha_levels):
        if abs(alev - alpha_eq) < 0.01:
            level_tags.append(f'α ≥ {alev:.2f} (equilibrium)')
        elif alev < alpha_eq:
            level_tags.append(f'α ≥ {alev:.2f} (below eq.)')
        else:
            level_tags.append(f'α ≥ {alev:.2f} (above eq.)')

    fig = plt.figure(figsize=(5.5 * n_cols, 5 * n_rows))

    for row, (alpha_level, times, tag) in enumerate(
            zip(alpha_levels, snap_times_sec_per_level, level_tags)):
        indices = [np.argmin(np.abs(t_sec - ts)) for ts in times]
        for col, idx in enumerate(indices):
            ax = fig.add_subplot(n_rows, n_cols,
                                 row * n_cols + col + 1, projection='3d')
            mask = alpha[:, :, :, idx] >= alpha_level

            if mask.any():
                xs = X3[mask]
                ys = Y3[mask]
                zs = Z3[mask]
                cs = T[:, :, :, idx][mask]
                sc = ax.scatter(xs, ys, zs, c=cs, cmap='coolwarm',
                                s=8, alpha=0.6,
                                vmin=Tmin_global, vmax=Tmax_global)
                fig.colorbar(sc, ax=ax, shrink=0.45, label='T [K]')
            else:
                ax.text2D(0.25, 0.5, 'α < level', transform=ax.transAxes)

            ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
            ax.set_title(f't = {t_sec[idx]:.0f} s\n{tag}', fontsize=9)

    fig.suptitle(f'3-D cure front  (α_eq ≈ {alpha_eq:.3f})', fontsize=14, y=1.01)
    fig.tight_layout()
    return fig


# ── 3-D cross-section triptych (xy, xz, yz through centre) ──────────
def plot_cross_sections_3D(x, y, z, t, U, alpha,
                           snap_times_T_sec=None,
                           snap_times_alpha_sec=None):
    """
    For selected times, show three orthogonal cross-sections of T and α
    through the domain centre: (x-y), (x-z), and (y-z).

    Parameters
    ----------
    snap_times_T_sec     : snapshot times [s] for the temperature figure
    snap_times_alpha_sec : snapshot times [s] for the α figure (defaults
                           to *snap_times_T_sec* when not given)

    Returns
    -------
    fig_T, fig_a
    """
    T     = U_to_T(U)
    t_sec = t * tdiff
    M     = U.shape[0]
    mid   = M // 2

    if snap_times_T_sec is None:
        snap_times_T_sec = np.linspace(0, t_sec[-1], 5).tolist()
    if snap_times_alpha_sec is None:
        snap_times_alpha_sec = snap_times_T_sec

    T_indices = [np.argmin(np.abs(t_sec - ts)) for ts in snap_times_T_sec]
    a_indices = [np.argmin(np.abs(t_sec - ts)) for ts in snap_times_alpha_sec]

    Tmin = T.min(); Tmax = T.max()

    XY_X, XY_Y = np.meshgrid(x, y, indexing='ij')
    XZ_X, XZ_Z = np.meshgrid(x, z, indexing='ij')
    YZ_Y, YZ_Z = np.meshgrid(y, z, indexing='ij')

    # ---- Temperature cross-sections ----
    n_cols_T = len(T_indices)
    fig_T, ax_T = plt.subplots(3, n_cols_T,
                               figsize=(3.8 * n_cols_T, 10),
                               constrained_layout=True, squeeze=False)
    for j, idx in enumerate(T_indices):
        # xy plane (z = mid)
        ax = ax_T[0, j]
        m1 = ax.pcolormesh(XY_X, XY_Y, T[:, :, mid, idx], shading='auto',
                           cmap='coolwarm', vmin=Tmin, vmax=Tmax)
        ax.set_xlabel('x'); ax.set_ylabel('y')
        ax.set_title(f't = {t_sec[idx]:.0f} s', fontsize=9)
        if j == 0:
            ax.set_ylabel('x-y plane\ny')
        fig_T.colorbar(m1, ax=ax, label='T [K]')

        # xz plane (y = mid)
        ax = ax_T[1, j]
        m1 = ax.pcolormesh(XZ_X, XZ_Z, T[:, mid, :, idx], shading='auto',
                           cmap='coolwarm', vmin=Tmin, vmax=Tmax)
        ax.set_xlabel('x'); ax.set_ylabel('z')
        if j == 0:
            ax.set_ylabel('x-z plane\nz')
        fig_T.colorbar(m1, ax=ax, label='T [K]')

        # yz plane (x = mid)
        ax = ax_T[2, j]
        m1 = ax.pcolormesh(YZ_Y, YZ_Z, T[mid, :, :, idx], shading='auto',
                           cmap='coolwarm', vmin=Tmin, vmax=Tmax)
        ax.set_xlabel('y'); ax.set_ylabel('z')
        if j == 0:
            ax.set_ylabel('y-z plane\nz')
        fig_T.colorbar(m1, ax=ax, label='T [K]')

    fig_T.suptitle('3-D Temperature – orthogonal cross-sections', fontsize=14)

    # ---- Alpha cross-sections ----
    n_cols_a = len(a_indices)
    fig_a, ax_a = plt.subplots(3, n_cols_a,
                               figsize=(3.8 * n_cols_a, 10),
                               constrained_layout=True, squeeze=False)
    for j, idx in enumerate(a_indices):
        ax = ax_a[0, j]
        m1 = ax.pcolormesh(XY_X, XY_Y, alpha[:, :, mid, idx], shading='auto',
                           cmap='Oranges', vmin=0, vmax=1)
        ax.set_xlabel('x'); ax.set_ylabel('y')
        ax.set_title(f't = {t_sec[idx]:.0f} s', fontsize=9)
        if j == 0:
            ax.set_ylabel('x-y plane\ny')
        fig_a.colorbar(m1, ax=ax, label='α')

        ax = ax_a[1, j]
        m1 = ax.pcolormesh(XZ_X, XZ_Z, alpha[:, mid, :, idx], shading='auto',
                           cmap='Oranges', vmin=0, vmax=1)
        ax.set_xlabel('x'); ax.set_ylabel('z')
        if j == 0:
            ax.set_ylabel('x-z plane\nz')
        fig_a.colorbar(m1, ax=ax, label='α')

        ax = ax_a[2, j]
        m1 = ax.pcolormesh(YZ_Y, YZ_Z, alpha[mid, :, :, idx], shading='auto',
                           cmap='Oranges', vmin=0, vmax=1)
        ax.set_xlabel('y'); ax.set_ylabel('z')
        if j == 0:
            ax.set_ylabel('y-z plane\nz')
        fig_a.colorbar(m1, ax=ax, label='α')

    fig_a.suptitle('3-D Curing (α) – orthogonal cross-sections', fontsize=14)

    return fig_T, fig_a


# ── Cure at centre & edge vs time (3-D) ──────────────────────────────
def plot_cure_center_3D(t, alpha):
    """
    Plot α at the 3-D domain centre, face centre, and corner over time.

    Returns
    -------
    fig, ax
    """
    t_sec = t * tdiff
    M     = alpha.shape[0]
    mid   = M // 2

    alpha_center = alpha[mid, mid, mid, :]
    alpha_face   = alpha[0, mid, mid, :]     # centre of one face
    alpha_corner = alpha[0, 0, 0, :]         # corner

    alpha_eq = _alpha_equilibrium()

    fig, ax = plt.subplots()
    ax.plot(t_sec, alpha_center, 'b-', linewidth=2, label='α at centre')
    ax.plot(t_sec, alpha_face,   'g-', linewidth=2, label='α at face centre')
    ax.plot(t_sec, alpha_corner, 'm-', linewidth=2, label='α at corner')
    ax.axhline(0.9, color='r', linestyle='--', linewidth=1, label='90 % cure')
    ax.axhline(alpha_eq, color='purple', linestyle=':', linewidth=1,
               label=f'equilibrium α ≈ {alpha_eq:.3f}')

    idx90 = np.where(alpha_center >= 0.9)[0]
    if len(idx90) > 0:
        t90 = t_sec[idx90[0]]
        ax.axvline(t90, color='gray', linestyle=':', linewidth=1)
        ax.annotate(f't₉₀ = {t90:.1f} s', xy=(t90, 0.9),
                    xytext=(t90 + 0.02 * t_sec[-1], 0.82),
                    arrowprops=dict(arrowstyle='->', color='gray'),
                    fontsize=10, color='gray')
        print(f"[3D] Time to 90 % cure at centre: {t90:.2f} s")
    else:
        print("[3D] 90 % cure NOT reached at the centre within simulation time.")

    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Degree of cure α')
    ax.set_title('3-D: Cure at centre, face, and corner')
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    return fig, ax


# ── Maximum temperature vs time (3-D) ────────────────────────────────
def plot_max_temp_3D(t, U, alpha=None):
    """
    Plot max(T) over the 3-D domain as a function of time.

    Returns
    -------
    fig, ax
    """
    t_sec = t * tdiff
    T     = U_to_T(U)
    T_max = T.max(axis=(0, 1, 2))

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
        print(f"[3D] Overheating!  Peak T = {T_max[peak_idx]:.2f} K "
              f"at t = {t_sec[peak_idx]:.2f} s")
    else:
        print("[3D] No overheating: max T never exceeds Tm.")

    y_lo = max(T_max.min(), Tm - 5)
    ax.set_ylim(y_lo, T_max.max() + 2)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Temperature [K]')
    ax.set_title('3-D: Maximum temperature inside the domain')
    ax.grid(True, alpha=0.3)

    if alpha is not None:
        M_  = alpha.shape[0]
        mid = M_ // 2
        ax2 = ax.twinx()
        ax2.plot(t_sec, alpha[mid, mid, mid, :], 'b-', linewidth=1.5,
                 label='α at centre')
        ax2.set_ylabel('Degree of cure α', color='b')
        ax2.set_ylim(0, 1.05)
        ax2.tick_params(axis='y', labelcolor='b')
        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, fontsize=8,
                  loc='center right')
    else:
        ax.legend(fontsize=8)

    return fig, ax


# ═══════════════════════════════════════════════════════════════════════
#  T > Tm  OVERHEATING  GIFS
# ═══════════════════════════════════════════════════════════════════════

# ── 2-D animated T > Tm ──────────────────────────────────────────────
def plot_T_above_Tm(x, y, t, U, max_frames=80):
    """
    Animate the temperature field, showing **only** regions where T > Tm.

    Frames are restricted to the time window during which at least one
    grid point exceeds Tm, so the GIF highlights *when* and *where*
    overheating occurs.

    Returns
    -------
    fig, ax, anim   (all *None* if T never exceeds Tm)
    """
    T     = U_to_T(U)
    t_sec = t * tdiff

    T_max_t = T.max(axis=(0, 1))          # max over space, per time step
    above   = T_max_t > Tm

    if not above.any():
        print("T never exceeds Tm — no overheating GIF to create.")
        return None, None, None

    above_idx = np.where(above)[0]
    start_idx = above_idx[0]
    end_idx   = above_idx[-1] + 1

    skip          = max(1, (end_idx - start_idx) // max_frames)
    frame_indices = np.arange(start_idx, end_idx, skip)

    X, Y = np.meshgrid(x, y, indexing='ij')
    vmin = Tm
    vmax = T[:, :, start_idx:end_idx].max()

    fig, ax = plt.subplots()
    data0 = T[:, :, frame_indices[0]].copy()
    data0[data0 <= Tm] = np.nan
    mesh = ax.pcolormesh(X, Y, data0, shading='auto',
                         cmap='hot', vmin=vmin, vmax=vmax)
    fig.colorbar(mesh, ax=ax, label='T [K]  (T > Tm only)')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    title = ax.set_title(f'T > Tm  |  t = {t_sec[frame_indices[0]]:.1f} s')

    def update(i):
        idx  = frame_indices[i]
        data = T[:, :, idx].copy()
        data[data <= Tm] = np.nan
        mesh.set_array(data.ravel())
        title.set_text(f'T > Tm  |  t = {t_sec[idx]:.1f} s')
        return mesh, title

    anim = FuncAnimation(fig, update, frames=len(frame_indices),
                         interval=100, blit=False)
    return fig, ax, anim


# ── 3-D animated T > Tm at z-slices ─────────────────────────────────
def plot_T_above_Tm_3D_slices(x, y, z, t, U,
                               z_fracs=(0.0, 0.25, 0.5),
                               max_frames=80):
    """
    Animate T(x, y, z_i, t) side-by-side for several z-slices, showing
    **only** regions where T > Tm.

    Returns
    -------
    fig, axes, anim   (all *None* if T never exceeds Tm)
    """
    T     = U_to_T(U)
    t_sec = t * tdiff
    M     = U.shape[0]

    z_indices = [min(int(f * (M - 1)), M - 1) for f in z_fracs]

    T_max_t = T.max(axis=(0, 1, 2))
    above   = T_max_t > Tm

    if not above.any():
        print("[3D] T never exceeds Tm — no overheating GIF to create.")
        return None, None, None

    above_idx = np.where(above)[0]
    start_idx = above_idx[0]
    end_idx   = above_idx[-1] + 1

    skip          = max(1, (end_idx - start_idx) // max_frames)
    frame_indices = np.arange(start_idx, end_idx, skip)

    X, Y = np.meshgrid(x, y, indexing='ij')
    vmin = Tm
    vmax = T[:, :, z_indices, start_idx:end_idx].max()

    n_slices = len(z_indices)
    fig, axes = plt.subplots(1, n_slices, figsize=(5 * n_slices, 4.5),
                             constrained_layout=True)
    if n_slices == 1:
        axes = [axes]

    meshes = []
    for ax, zi, zf in zip(axes, z_indices, z_fracs):
        data0 = T[:, :, zi, frame_indices[0]].copy()
        data0[data0 <= Tm] = np.nan
        m1 = ax.pcolormesh(X, Y, data0, shading='auto',
                           cmap='hot', vmin=vmin, vmax=vmax)
        fig.colorbar(m1, ax=ax, label='T [K]')
        ax.set_xlabel('x'); ax.set_ylabel('y')
        ax.set_title(f'z = {z[zi]:.2f}')
        meshes.append(m1)

    suptitle = fig.suptitle(
        f'T > Tm only  |  t = {t_sec[frame_indices[0]]:.1f} s', fontsize=13)

    def update(i):
        idx = frame_indices[i]
        for m1, zi in zip(meshes, z_indices):
            data = T[:, :, zi, idx].copy()
            data[data <= Tm] = np.nan
            m1.set_array(data.ravel())
        suptitle.set_text(f'T > Tm only  |  t = {t_sec[idx]:.1f} s')
        return meshes + [suptitle]

    anim = FuncAnimation(fig, update, frames=len(frame_indices),
                         interval=100, blit=False)
    return fig, axes, anim




def plot_error_convergence(h, errs, degree,coarse=False):
    """
    Log-log plot of max-norm error vs step size with fitted convergence order.

    Returns
    -------
    fig, ax
    """
    fig, ax = plt.subplots()
    ax.plot(h, errs, 'o-', label=f'Approximate order = {degree:.2f}')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Step size', fontsize=13)
    ax.set_ylabel(r'$\|E\|_{\infty}$', fontsize=13)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)
    if coarse:
        ax.set_title(r'Error convergence by grid refinement $N = 5M$')
    else:
        ax.set_title(r'Error convergence by grid refinement $N = 10M$')
    return fig, ax

