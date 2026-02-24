## Task 1

We have the following model:
$$
\begin{align}
\rho c_p \frac{\partial T}{\partial t} &= k\Delta T + \rho H_r \frac{\partial \alpha}{\partial t} \quad \text{on} \; \Omega, \\
-k\frac{\partial T}{\partial n} &= h(T-T_m) \quad\quad\quad \text{on} \;\partial\Omega.
\end{align}
$$

To rewrite the model in a dimensionless form, we use the scaling suggested in the task:
$$
\begin{align*}
u = \frac{T-T_m}{T_m-T_0}, \quad x^* = x/L, \quad y^* = y/L \\ t^* =\frac{t}{t_\text{diff}}, \quad t_\text{diff} = \frac{\rho c_pL^2}{k}.
\end{align*}
$$

We now need to find
$$
\frac{\partial T}{\partial t}, \;\Delta T , \;\frac{\partial \alpha}{\partial t} \;\text{and}\; \frac{\partial T}{\partial n},
$$
expressed with
$$
\frac{\partial u}{\partial t^*}, \;\Delta^*u = \frac{\partial^2 u}{\partial (x^*)^2} + \frac{\partial^2 u}{\partial (y^*)^2},\, \;\frac{\partial \alpha}{\partial t^*} \;\text{and}\; \frac{\partial u}{\partial n^*},
$$
and substitute them into the model. We use 
$$
\begin{align*}
\frac{\partial u}{\partial T} &= \frac{1}{T_m - T_0}, \\
\frac{\partial t^*}{\partial t} &= \frac{1}{t_\text{diff}}, \\
\frac{\partial x^*}{\partial x} = \frac{\partial y^*}{\partial y} &= \frac{1}{L},
\end{align*}
$$
and get the following relationships: 
$$
\begin{align*}
\frac{\partial T}{\partial t} &= \frac{T_m-T_0}{t_\text{diff}}u_{t^*},\quad u_{t^*} = \frac{\partial u}{\partial t^*} \\
\Delta T &= \frac{T_m-T_0}{L^2}\Delta^*u, \\
\frac{\partial \alpha}{\partial t} &= \frac{1}{t_\text{diff}}\alpha_{t^*}, \quad\quad\quad\alpha_{t^*} = \frac{\partial\alpha}{\partial t^*}, \\
\frac{\partial T}{\partial n} &= \frac{T_m - T_0}{L}u_{n^*}, \quad u_{n^*} = \frac{\partial u}{\partial n^*}.
\end{align*}
$$

We plug this into the model on $\Omega$, eq. $(1)$, and get:
$$
\begin{align*}
\rho c_p\frac{\partial T}{\partial t} &= k\Delta T + \rho H_r \frac{\partial \alpha}{\partial t} \\
\rho c_p \frac{T_m-T_0}{t_\text{diff}}u_{t^*} &= k\frac{T_m-T_0}{L^2}\Delta^*u + \rho H_r \frac{1}{t_\text{diff}}\alpha_{t^*}\\ 
\rho c_p \frac{T_m-T_0}{t_\text{diff}}u_{t^*} &= \rho c_p\frac{T_m-T_0}{t_\text{diff}}\Delta^*u + \rho H_r \frac{1}{t_\text{diff}}\alpha_{t^*}\\
u_{t^*} &= \Delta^*u + \frac{H_r}{c_p(T_m-T_0)} \alpha_{t^*}\\  
u_{t^*} &= \Delta^*u + c \alpha_{t^*}, \quad c= \frac{H_r}{c_p(T_m-T_0)}. 
\end{align*}
$$ 
This equation is dimensionless, which is what we wanted. We do the same for the model on $\partial\Omega$, eq. $(2)$, and get:
$$
\begin{align*}
-k\frac{\partial T}{\partial n} &= h(T-T_m) \\
-k\frac{T_m - T_0}{L}u_{n^*} &= h(T-T_m) \\
u_{n^*} &= -\frac{hL}{k}u \\
u_{n^*} + \gamma u &= 0, \quad \gamma=\frac{hL}{k}.
\end{align*}
$$
Now we have the model in dimensionless form.