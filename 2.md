We use the following approximation for dicretizations with known error term: 
$$
\begin{aligned}
\frac{1}{k}\nabla_k U_p =  \frac{U_p^n - U_p^{n-1}}{k} = u_t -\frac{k}{2}u_{xx}(\xi_m)\\
\frac{1}{h^2}(\delta_x^2 + \delta_y^2) = \Delta_* u +\frac{h^2}{12}(u_{4x}(\xi_m)+u_{4y}(\xi_m))
\end{aligned}
$$

Inserting the discretizations in to our dimensionless equation yields the following: 

$$
\begin{aligned} 
\frac{U_p^{n} - U_p^{n-1}}{k} = \frac{1}{h^2}(-4U_P^n+U_W^n+U_E^n+U_S^n+U_N^n) +kc\alpha_{t^*} \\
U_p^{n} = U_p^{n-1} + \frac{k}{h^2}(-4U_P^n+U_W^n+U_E^n+U_S^n+U_N^n) +kc \alpha_{t^*} \\ 
(4+\frac{h^2}{k})U_p^n -U_W^n -U_N^n-U_S^n-U_E^n = \frac{h^2}{k}U_p^{n-1} + h^2c\alpha_{t^*}\\ 
\end{aligned} 
$$



Which gives us the following linear system for the interior points: 
$$
\begin{aligned}
A_hU_p^n = b\\
A = \text{tridiag}(-I_{M-1},B,-I_{M-1}),\quad B = \text{tridiag}(-1,4+\frac{h^2}{k},-1)\\
b = \frac{h^2}{k}U_p^{n-1} + h^2c \alpha_{t^*}\\ 
\end{aligned} 
$$



$$
\begin{align}

\end{align}
$$

For stability and consistency analysis we ignore the source term.
Truncation error is given by: 
$$
\begin{aligned}
\tau = (L-L_h)u_p = (\partial t-\frac{1}{k}\nabla_k +\frac{1}{h^2}(\delta_x^2 + \delta_y^2)-\Delta)u_p\\ 
\tau = \frac{k}{2}u_{tt}(\lambda_n) + \frac{h^2}{12}(u_{4x}(\xi_m)+u_{4y}(\eta_m))\\
\tau \le O(k + h^2)\quad, \lambda_n \in [t_n,t_{n+1}], \xi_m \in[x_m,x_{m+1}], \eta_m \in [y_m,y_{m+1}] \
\end{aligned}
$$
Thus we have shown that our method is consistent and of first order in time and second order in space. 

For the Von Neumann stability analysis we put $U_p^n = \xi^n e^{ib0}$ for our numerical scheme and obtain: 
$$\xi^n = \xi^{n-1} + \frac{k}{h^2}\xi^n(-4 + e^{ib_xh} +  e^{-ib_xh} + e^{ib_yh} +  e^{-ib_yh}) \xi= 1 + \frac{2k}{h^2}(-2\sin^2(\frac{b_xh}{2})-2\sin^2(\frac{b_yh}{2}))\\ 
\quad \\ \xi = \frac{1}{1 + \frac{4k}{h^2}(2\sin^2(\frac{b_xh}{2})+2\sin^2(\frac{b_yh}{2}))} \\ \quad\\ |\xi| \le 1 \quad  \forall \beta.$$

And thus our method is Von Neumann stable, which is expected for backward difference discretization for the heat equation in time.