"""
I/O and caching utilities.
"""

import os
import numpy as np

from src.solvers import solve_heat, solve_vulcanization

# Default output root (relative to repo root)
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "output")


def _ensure_dirs():
    """Create output sub-directories if they don't exist."""
    for sub in ("solutions", "plots"):
        os.makedirs(os.path.join(OUTPUT_DIR, sub), exist_ok=True)


def _cache_path(problem, M, N, T_end):
    return os.path.join(OUTPUT_DIR, "solutions",
                        f"{problem}_M{M}_N{N}_T{T_end}.npz")


# ── Save / load raw arrays ────────────────────────────────────────────
def save_solution(path, x, y, t, U, alpha=None):
    """Save solver output to a compressed .npz file."""
    data = dict(x=x, y=y, t=t, U=U)
    if alpha is not None:
        data["alpha"] = alpha
    np.savez_compressed(path, **data)
    print(f"Solution saved to {path}")


def load_solution(path):
    """
    Load a solution saved with :func:`save_solution`.

    Returns ``(x, y, t, U)`` or ``(x, y, t, U, alpha)``.
    """
    data = np.load(path)
    x, y, t, U = data["x"], data["y"], data["t"], data["U"]
    if "alpha" in data:
        return x, y, t, U, data["alpha"]
    return x, y, t, U


# ── High-level: solve-or-load ─────────────────────────────────────────
def get_solution(problem, M=51, N=101, T_end=None, force=False):
    """
    Return a cached solution or solve and cache it.

    Parameters
    ----------
    problem : ``"Heat"`` or ``"Vulcanization"``
    force   : if True, always re-solve
    """
    _ensure_dirs()

    if T_end is None:
        T_end = 1 if problem == "Heat" else 20

    path = _cache_path(problem, M, N, T_end)

    if not force and os.path.exists(path):
        print(f"Loading cached solution from {path}")
        return load_solution(path)

    print(f"Solving {problem} (M={M}, N={N}, T_end={T_end}) …")
    if problem == "Heat":
        sol = solve_heat(M=M, N=N, T_end=T_end)
        save_solution(path, *sol)
        return sol
    elif problem == "Vulcanization":
        sol = solve_vulcanization(M=M, N=N, T_end=T_end)
        save_solution(path, *sol[:4], alpha=sol[4])
        return sol
    else:
        raise ValueError(f"Unknown problem: {problem!r}")


def plots_dir():
    """Return the absolute path to the plots output directory."""
    _ensure_dirs()
    return os.path.join(OUTPUT_DIR, "plots")
