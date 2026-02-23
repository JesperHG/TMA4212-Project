import numpy as np
import matplotlib.pyplot as plt
from Heat import Heat

from MODEL_PARAMS import k, h, rho, c_p, H_r, A, E, R, m, n, T_s, L, T0, Tm


f = lambda x, y: 0
g = lambda x, y: -np.exp(-h*L/k)

pde = Heat(f,g)

solves = []
errs = []
M = np.array([10,20,40,80,160])

for i in M:
    x, y, t, U = pde.solve(M=int(i),N=50)
    
    solves.append(U)
    if len(solves) > 1:
        err = solves[-1][::2,::2] - solves[-2]
        errs.append(np.max(np.abs(err)))

errs = np.array(errs)

h = 1/M[1:]

degree = np.polyfit(np.log(h),np.log(errs),1)[0]
print(degree)


plt.loglog(h,errs)
plt.show()






