import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from MODEL_PARAMS import k, h, rho, c_p, H_r, A, E, R, m, n, T_s, L, T0, Tm
from Heat import Heat
from Vulcanization import Vulcanization
from plotting import plotU

pde = Heat()
pde2 = Vulcanization(0,0)

x, y, t, U, alpha = pde2.solve() 

# Create animated plot
fig, ax, anim = plotU(x, y, t, U)
plt.show()  

fig, ax, anim = plotU(x, y, t, alpha)
plt.show()  



