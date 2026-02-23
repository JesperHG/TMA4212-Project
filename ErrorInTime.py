
import numpy as np
import matplotlib.pyplot as plt
from Heat import Heat

from MODEL_PARAMS import k, h, rho, c_p, H_r, A, E, R, m, n, T_s, L, T0, Tm
f = lambda x, y: 0
g = lambda x, y: -np.exp(-h*L/k)

pde = Heat(f,g)

solves = []
errs = []

N = np.array([4,8,16,32,64])*100

for j in N:
    x, y, t, U = pde.solve(M=10,N=int(j))
    
    solves.append(U)
    if len(solves) > 1:
        err = solves[-1][::1,::1,::2] - solves[-2]
        errs.append(np.max(np.abs(err)))

errs = np.array(errs)

k = 1/N[1:]

degree = np.polyfit(np.log(k),np.log(errs),1)[0]
print(degree)


plt.loglog(k,errs)
plt.show()