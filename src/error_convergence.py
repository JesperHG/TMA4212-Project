"""
Error convergence analysis for the heat equation.
"""

import numpy as np
from src.utils import get_solution


def error_test(steps, param="pos"):
    """
    Convergence test for the heat equation via grid refinement.

    Computes the max-norm error between a fine-grid reference solution
    and coarser solutions for each value in *steps*. A log-log linear
    fit gives the approximate convergence order.

    Parameters
    ----------
    steps : array_like of int
        Grid sizes to test. Interpreted as spatial grid counts M when
        ``param="pos"``, or temporal grid counts N when ``param="time"``.
    param : {"pos", "time"}, optional
        Which discretisation parameter to refine.
        ``"pos"``  – refine in space (M varies, N fixed).
        ``"time"`` – refine in time  (N varies, M fixed).

    Returns
    -------
    errs : ndarray, shape (len(steps),)
        Max-norm errors for each step size.
    h : ndarray, shape (len(steps),)
        Corresponding step sizes (1 / steps).
    degree : float
        Estimated convergence order from the log-log fit.
    """
    n_steps = len(steps)
    errs = np.zeros(n_steps)

    if param == "pos": 
        N_fixed = 50
        U_ref = get_solution("Heat", M=1000, N=N_fixed)[-1]
        U_ref_spaceavg = np.mean(U_ref, axis=(0, 1)) # Avg over first two axes (spatial)

        for i in range(n_steps):
            U = get_solution("Heat", M=steps[i], N=N_fixed)[-1]
            U_spaceavg = np.mean(U, axis=(0, 1)) # Avg over first two axes (spatial)

            errs[i] = np.max(np.abs(U_ref_spaceavg - U_spaceavg))

    elif param == "time": 
        M_fixed=20
        U_ref = get_solution("Heat", M=M_fixed, N=10000)[-1]
        U_ref_timeavg = np.mean(U_ref, axis=-1) # Avg over last axis (time)

        for i in range(n_steps):
            U = get_solution("Heat", M=M_fixed, N=steps[i])[-1]
            U_timeavg = np.mean(U, axis=-1) # Avg over last axis (time)

            errs[i] = np.max(np.abs(U_ref_timeavg - U_timeavg))

    h = 1 / steps
    degree = np.polyfit(np.log(h), np.log(errs), 1)[0]
    return errs, h, degree

