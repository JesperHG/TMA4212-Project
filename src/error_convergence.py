"""
Error convergence analysis for the heat equation.
"""

import numpy as np
from src.solvers import solve_heat


def error_test(M_steps, N_steps, param="pos"):
    """
    Run a sequence of heat-equation solves with increasing resolution and
    compute the max-norm error between successive refinements.

    Parameters
    ----------
    M_steps : array-like of spatial grid sizes (each doubling)
    N_steps : array-like of temporal grid sizes (each doubling)
    param   : ``"pos"`` for spatial refinement, ``"time"`` for temporal

    Returns
    -------
    errs   : 1-D array of max-norm errors
    h      : 1-D array of step sizes
    degree : fitted convergence order (slope in log-log)
    """
    solves = []
    errs = []

    if param == "pos":
        for i in range(len(M_steps)):
            x, y, t, U = solve_heat(int(M_steps[i]), int(N_steps[i]))
            solves.append(U)
            if len(solves) > 1:
                err = solves[-1][::2, ::2, -1] - solves[-2][:, :, -1]
                errs.append(np.max(np.abs(err)))

        h = 1 / M_steps[1:]
        errs = np.array(errs)

    elif param == "time":
        for i in range(len(M_steps)):
            x, y, t, U = solve_heat(int(M_steps[i]), int(N_steps[i]))
            solves.append(U)
            if len(solves) > 1:
                err = solves[-1][:, :, -1] - solves[-2][:, :, -1]
                errs.append(np.max(np.abs(err)))

        h = 1 / N_steps[1:]
        errs = np.array(errs)

    degree = np.polyfit(np.log(h), np.log(errs), 1)[0]
    return errs, h, degree

