"""
Main script – solve and generate all plots for the project.

Run from the project root:
    python main.py
"""

import os
import matplotlib.pyplot as plt

from src.plotting import (plot_temperature, plot_alpha,
                          plot_T_snapshots, plot_alpha_snapshots,
                          plot_cure_center, plot_max_temp,
                          plot_temperature_3D_slices, plot_alpha_3D_slices,
                          plot_T_snapshots_3D, plot_alpha_snapshots_3D,
                          plot_centerline_z, plot_isosurface_3D,
                          plot_cross_sections_3D,
                          plot_cure_center_3D, plot_max_temp_3D,
                          plot_T_above_Tm, plot_T_above_Tm_3D_slices)
from src.utils import get_solution, plots_dir


# ── Helpers ────────────────────────────────────────────────────────────
def _save(fig, subfolder, name, dpi=200):
    path = os.path.join(plots_dir(subfolder), name)
    fig.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {path}")


def _save_anim(anim, subfolder, name, fps=10):
    path = os.path.join(plots_dir(subfolder), name)
    anim.save(path, writer='pillow', fps=fps)
    plt.close(anim._fig)
    print(f"  Saved {path}")


# ── Main ──────────────────────────────────────────────────────────────
def main():
    # === Pure heat equation ================================================
    xh, yh, th, Uh = get_solution("Heat", M=51, N=2001, T_end=1)

    fig, ax, anim = plot_temperature(xh, yh, th, Uh, max_frames=60, t_end_sec=1000)
    _save_anim(anim, "gifs", "heat_T.gif", fps=8)

    # === Vulcanization (heat + reaction) ===================================
    x, y, t, U, alpha = get_solution("Vulcanization", M=51, N=2001, T_end=20)

    # -- GIFs ---------------------------------------------------------------
    fig, ax, anim = plot_temperature(x, y, t, U, max_frames=100, t_end_sec=1000)
    _save_anim(anim, "gifs", "vulc2d_T.gif", fps=8)

    fig, ax, anim = plot_alpha(x, y, t, alpha, max_frames=40)
    _save_anim(anim, "gifs", "vulc2d_alpha.gif", fps=8)

    # T > Tm overheating GIF (2-D)
    result = plot_T_above_Tm(x, y, t, U, max_frames=80)
    if result[2] is not None:
        _save_anim(result[2], "gifs", "vulc2d_T_above_Tm.gif", fps=8)

    # -- Snapshots ----------------------------------------------------------
    fig, axes = plot_T_snapshots(
        x, y,
        [th, t], [Uh, U],
        labels=['Heat', 'Vulcanization'],
        snap_times_sec=[0, 100, 200, 500, 1000])
    _save(fig, "snapshots", "heat_vs_vulc_T.png")

    fig, axes = plot_alpha_snapshots(
        x, y, t, alpha,
        snap_times_sec=[0, 1000, 20000, 40000, 60000, 80000])
    _save(fig, "snapshots", "vulc2d_alpha.png")

    # -- Analysis -----------------------------------------------------------
    fig, ax = plot_cure_center(t, alpha)
    _save(fig, "analysis", "vulc2d_cure_center.png")

    fig, ax = plot_max_temp(t, U, alpha=alpha)
    _save(fig, "analysis", "vulc2d_max_T.png")

    # ═══════════════════════════════════════════════════════════════════════
    #  3-D Vulcanization
    # ═══════════════════════════════════════════════════════════════════════
    print("\n=== 3-D Vulcanization ===")
    x3, y3, z3, t3, U3, alpha3 = get_solution(
        "Vulcanization3D", M=15, N=2001, T_end=20)

    z_slices = (0.0, 0.25, 0.5)  # boundary, quarter, centre

    # -- GIFs ---------------------------------------------------------------
    fig, axes, anim = plot_temperature_3D_slices(
        x3, y3, z3, t3, U3, z_fracs=z_slices,
        max_frames=80, t_end_sec=1000)
    _save_anim(anim, "gifs", "vulc3d_T_slices.gif", fps=8)

    fig, axes, anim = plot_alpha_3D_slices(
        x3, y3, z3, t3, alpha3, z_fracs=z_slices, max_frames=60)
    _save_anim(anim, "gifs", "vulc3d_alpha_slices.gif", fps=8)

    # T > Tm overheating GIF (3-D z-slices)
    result = plot_T_above_Tm_3D_slices(
        x3, y3, z3, t3, U3, z_fracs=z_slices, max_frames=80)
    if result[2] is not None:
        _save_anim(result[2], "gifs", "vulc3d_T_above_Tm_slices.gif", fps=8)

    # -- Snapshots ----------------------------------------------------------
    fig, axes = plot_T_snapshots_3D(
        x3, y3, z3, t3, U3, z_fracs=z_slices,
        snap_times_sec=[0, 100, 200, 500, 1000])
    _save(fig, "snapshots", "vulc3d_T_slices.png")

    fig, axes = plot_alpha_snapshots_3D(
        x3, y3, z3, t3, alpha3, z_fracs=z_slices,
        snap_times_sec=[0, 1000, 20000, 40000, 60000, 80000])
    _save(fig, "snapshots", "vulc3d_alpha_slices.png")

    # Orthogonal cross-sections (separate T and α times)
    fig_T, fig_a = plot_cross_sections_3D(
        x3, y3, z3, t3, U3, alpha3,
        snap_times_T_sec=[0, 200, 500, 1000],
        snap_times_alpha_sec=[0, 20000, 40000, 80000])
    _save(fig_T, "snapshots", "vulc3d_cross_T.png")
    _save(fig_a, "snapshots", "vulc3d_cross_alpha.png")

    # 3-D isosurface scatter (3×3: alpha levels × time)
    fig = plot_isosurface_3D(x3, y3, z3, t3, U3, alpha3)
    _save(fig, "snapshots", "vulc3d_isosurface.png")

    # -- Analysis -----------------------------------------------------------
    # Centre-line profiles along z (separate T and α times)
    fig, axes = plot_centerline_z(
        z3, t3, U3, alpha=alpha3,
        snap_times_T_sec=[0, 50, 100, 200, 500, 1000],
        snap_times_alpha_sec=[0, 5000, 20000, 40000, 60000, 80000])
    _save(fig, "analysis", "vulc3d_centerline_z.png")

    fig, ax = plot_cure_center_3D(t3, alpha3)
    _save(fig, "analysis", "vulc3d_cure_center.png")

    fig, ax = plot_max_temp_3D(t3, U3, alpha=alpha3)
    _save(fig, "analysis", "vulc3d_max_T.png")

    print(f"\nAll plots saved to {plots_dir()}/")


if __name__ == "__main__":
    main()
