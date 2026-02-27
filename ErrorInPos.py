import numpy as np
import matplotlib.pyplot as plt
from src.solvers import solve_heat




solves = []
errs = []
M = np.array([100,200,400,800])

for i in M:
    x, y, t, U = solve_heat(128,N=int(i))

    solves.append(U)
    if len(solves) > 1:
        err = solves[-1][:,:,-1] - solves[-2][:,:,-1]
        errs.append(np.max(np.abs(err)))

errs = np.array(errs)

h = 1/M[1:]

degree = np.polyfit(np.log(h),np.log(errs),1)[0]
print(degree)


plt.loglog(h,errs)
plt.show()






