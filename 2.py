import numpy as np
import matplotlib.pyplot as plt

from MODEL_PARAMS import k, h, rho, c_p, H_r, A, E, R, m, n, T_s
from Heat import Heat

L = 0.02 # [m] Domain size
T0 = 300 # [K] Initial temperature
Tm = 443 # [K] Mold temperature

f = lambda x, y: 0
g = lambda x, y: -np.exp(-k*L/h)

pde = Heat(f, g)
x, y, U = pde.solve(4)