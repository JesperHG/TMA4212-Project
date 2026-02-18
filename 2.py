import numpy as np
import matplotlib.pyplot as plt

L = 0.02 # [m] Domain size
T0 = 300 # [K] Initial temperature
Tm = 443 # [K] Mold temperature

Omega = np.linspace((0,0), (L,L), 100) # Square domain
