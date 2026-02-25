"""
Task 2 – Vulcanisation simulation: solve, generate and save all plots.

Run from the project root:
    python -m tasks.task2
"""

import os
import matplotlib.pyplot as plt

from src.plotting import (plot_temperature, plot_alpha, plot_snapshots,
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
    x, y, t, U = get_solution("Heat", M=51, N=101, T_end=1)

    fig, ax, anim = plot_temperature(x, y, t, U, max_frames=100, t_end_sec=1000)
    _save_anim(anim, "heat_T.gif", fps=8)

    # --- Vulcanization (heat + reaction) ---
    x, y, t, U, alpha = get_solution("Vulcanization", M=51, N=1001, T_end=20)

    # Animated temperature field (only first ~1000 s)
    fig, ax, anim = plot_temperature(x, y, t, U, max_frames=100, t_end_sec=1000)
    _save_anim(anim, "vulc_T.gif", fps=8)

    # Animated curing field (sub-sampled)
    fig, ax, anim = plot_alpha(x, y, t, alpha, max_frames=40)
    _save_anim(anim, "vulc_alpha.gif", fps=8)

    # Snapshots at selected physical times (seconds)
    fig, axes = plot_snapshots(x, y, t, U, alpha,
                               snap_times_sec=[100, 200, 500, 5000, 20000, 80000])
    _save(fig, "vulc_snapshots.png")

    # Cure at centre vs time (with 90 % marker)
    fig, ax = plot_cure_center(t, alpha)
    _save(fig, "vulc_cure_center.png")

    # Max temperature vs time (with mold-temp marker)
    fig, ax = plot_max_temp(t, U)
    _save(fig, "vulc_max_temp.png")

    print(f"\nAll plots saved to {plots_dir()}/")


if __name__ == "__main__":
    main()
