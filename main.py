"""
Main script – solve and generate all plots for the project.

Run from the project root:
    python main.py
"""

import os
import matplotlib.pyplot as plt

from src.plotting import (plot_temperature, plot_alpha,
                          plot_T_snapshots, plot_alpha_snapshots,
                          plot_cure_center, plot_max_temp)
from src.utils import get_solution, plots_dir


# ── Helpers ────────────────────────────────────────────────────────────
def _save(fig, name, dpi=200):
    path = os.path.join(plots_dir(), name)
    fig.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {path}")


def _save_anim(anim, name, fps=10):
    path = os.path.join(plots_dir(), name)
    anim.save(path, writer='pillow', fps=fps)
    plt.close(anim._fig)
    print(f"  Saved {path}")


# ── Main ──────────────────────────────────────────────────────────────
def main():
    # --- Pure heat equation ---
    xh, yh, th, Uh = get_solution("Heat", M=51, N=1001, T_end=1)

    fig, ax, anim = plot_temperature(xh, yh, th, Uh, max_frames=60, t_end_sec=1000)
    _save_anim(anim, "heat_T.gif", fps=8)

    # --- Vulcanization (heat + reaction) ---
    x, y, t, U, alpha = get_solution("Vulcanization", M=51, N=1001, T_end=20)

    # Animated temperature field (only first ~1000 s)
    fig, ax, anim = plot_temperature(x, y, t, U, max_frames=100, t_end_sec=1000)
    _save_anim(anim, "vulc_T.gif", fps=8)

    # Animated curing field (sub-sampled)
    fig, ax, anim = plot_alpha(x, y, t, alpha, max_frames=40)
    _save_anim(anim, "vulc_alpha.gif", fps=8)

    # Temperature snapshots: heat vs vulcanization, early times
    fig, axes = plot_T_snapshots(
        x, y,
        [th, t], [Uh, U],
        labels=['Heat', 'Vulcanization'],
        snap_times_sec=[0, 100, 200, 500, 1000])
    _save(fig, "snapshots_T.png")

    # Curing snapshots: spread over full simulation
    fig, axes = plot_alpha_snapshots(
        x, y, t, alpha,
        snap_times_sec=[0, 1000, 20000, 40000, 60000, 80000])
    _save(fig, "snapshots_alpha.png")

    # Cure at centre & edge vs time (with 90 % marker + equilibrium)
    fig, ax = plot_cure_center(t, alpha)
    _save(fig, "vulc_cure_center.png")

    # Max temperature vs time (with cure on twin axis)
    fig, ax = plot_max_temp(t, U, alpha=alpha)
    _save(fig, "vulc_max_temp.png")

    print(f"\nAll plots saved to {plots_dir()}/")


if __name__ == "__main__":
    main()
