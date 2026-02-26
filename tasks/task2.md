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

Now to derive the linear system for the boundary points. We have collected the following equations:

$$
\begin{align*} 
4+\frac{h^2}{k}U_p^n -U_W^n -U_N^n-U_S^n-U_E^n &= \frac {h^2}{k}U_p^{n-1} + h^2c\alpha_{t^*} \quad \text{Interior}\\
\frac{1}{2h}(u_E - u_W) &= -\gamma u  \quad \text{Left and Right boundary}\\  

\frac{1}{2h}(u_N - u_S) &= -\gamma u  \quad \text{Upper and Bottom boundary}\\ 
\end{align*}
$$

Now solving for the different ghost points that occurs at the different boundaries gives us the following equations for the different boundaries:


$$
\begin{align*}
(4 + \frac{h^2}{k}+2h\gamma)U_p^n -2U_E^n -U_S^n -U_N^n  &= \frac{h^2}{k}U_p^{n-1} \quad \text{Left boundary}\\  
(4 + \frac{h^2}{k}+2h\gamma)U_p^n -2U_W^n -U_S^n -U_N^n  &= \frac{h^2}{k}U_p^{n-1} \quad \text{Right boundary}\\ 
(4 + \frac{h^2}{k}+2h\gamma)U_p^n -2U_N^n -U_E^n -U_W^n  &= \frac{h^2}{k}U_p^{n-1} \quad \text{Bottom boundary}\\
(4 + \frac{h^2}{k}+2h\gamma)U_p^n -2U_S^n -U_E^n -U_W^n  &= \frac{h^2}{k}U_p^{n-1} \quad \text{Upper boundary}
\end{align*}
$$ 

The only thing remaining is to take into account is the corner cases. For the corner cases we have to use both the horizontal and vertical boundary equations to remedy the two occuring ghost points. This yields the following  set of equations:
$$
\begin{align*}
(4 + \frac{h^2}{k}+4h\gamma)U_p^n -2U_E^n -2U_S^n   &= \frac{h^2}{k}U_p^{n-1} \quad \text{Upper left boundary}\\  
(4 + \frac{h^2}{k}+4h\gamma)U_p^n -2U_W^n -2U_S^n   &= \frac{h^2}{k}U_p^{n-1} \quad \text{Upper right boundary}\\ 
(4 + \frac{h^2}{k}+4h\gamma)U_p^n -2U_N^n -2U_E^n  &= \frac{h^2}{k}U_p^{n-1} \quad \text{Lower left boundary}\\
(4 + \frac{h^2}{k}+4h\gamma)U_p^n -2U_N^n -2U_W^n  &= \frac{h^2}{k}U_p^{n-1} \quad \text{Lower right boundary}
\end{align*}
$$  

Now to derive the linear system:  

$$
A_hU_p^n = b = \frac{h^2}{k}U_p^{n-1} + h^2c\alpha_{t^*} \\  
$$


First let us consider the system of equations for the points on the upper and lower boundary. The set of equations for the $M$ number of points at each boundary can be expressed as: 

$$
\begin{align*}
D &= 
\begin{bmatrix}
d^* & -2 & 0 & \cdots & 0 \\
-1 & d & -1 & \cdots & 0 \\
0 & -1 & d  & \ddots & \vdots \\
\vdots & \vdots & \ddots & \ddots & -1 \\
0 & 0 & \cdots & -2 & d^* 
\end{bmatrix}
\\
d^* &= 4 +4h\gamma + \frac{h^2}{k}\\
d &= 4 +2h\gamma + \frac{h^2}{k}
\end{align*}
$$   

Next to consider is the system of equations for the interior and the left and right boundary points which can be express as:

$$
\begin{align*}
D_{inner} &= 
\begin{bmatrix}
d & -2 & 0 & \cdots & 0 \\
-1 & c & -1 & \cdots & 0 \\
0 & -1 & c & \ddots & \vdots \\
\vdots & \vdots & \ddots & \ddots & -1 \\
0 & 0 & \cdots & -2 & d 
\end{bmatrix}
\\
c &= 4 + \frac{h^2}{k}\\
d &= 4 +2h\gamma + \frac{h^2}{k}
\end{align*}
$$    

Using this  we can construct the linear system for the whole grid as follows: 

$$
\begin{align*}
A &= 
\begin{bmatrix}
D & -2I & 0 & \cdots & 0 \\
-I & D_{inner} & -I & \cdots & 0 \\
0 & -I & D_{inner} & \ddots & \vdots \\
\vdots & \vdots & \ddots & \ddots & -I \\
0 & 0 & \cdots & -2I & D
\end{bmatrix}
\\ 
AU_p^{n} &= U_p^{n-1} + h^2c \alpha_{t^*}
\end{align*}
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
Thus we have shown that our method is consistent for the interior point and of first order in time and second order in space. 

For the Von Neumann stability analysis we put $U_p^n = \xi^n e^{i\beta0}$ for our numerical scheme. First we peform stability analysis on the interior equations:  
$$
U_p^{n} = U_p^{n-1} + \frac{k}{h^2}(-4U_P^n+U_W^n+U_E^n+U_S^n+U_N^n) \\ 
\xi^n = \xi^{n-1} + \frac{k}{h^2}\xi^n(-4 + e^{i\beta_xh} +  e^{-i\beta_xh} + e^{i\beta_yh} +  e^{-i\beta_yh}) \xi= 1 + \frac{2k}{h^2}(-2\sin^2(\frac{\beta_xh}{2})-2\sin^2(\frac{\beta_yh}{2}))\\ 
\quad \\ \xi = \frac{1}{1 + \frac{4k}{h^2}(2\sin^2(\frac{b_xh}{2})+2\sin^2(\frac{b_yh}{2}))} \\ \quad\\ |\xi| \le 1 \quad  \forall \beta_x,\beta_y.$$

Next is the boundary equations, we start by considering the right boundary equation:
$$
(4 + \frac{h^2}{k}+2h\gamma)U_p^n -2U_W^n -U_S^n -U_N^n = \frac{h^2}{k}U_p^{n-1}\\
\quad\\
(4 + \frac{h^2}{k} + 2h\gamma -2e^{-i\beta_xh} -e^{-i\beta_yh}-e^{i\beta_yh})\xi = \frac{h^2}{k}\\  
\quad\\
(2 + \frac{h^2}{k} + 2h\gamma -2e^{-i\beta_xh} + 4\sin^2(\beta_y h))\xi = \frac{h^2}{k}\\ 
\xi = \frac{1}{1+\frac{2k\gamma}{h} + 2\frac{k}{h^2} - \frac{2k}{h^2} e^{-i\beta_xh}+\frac{4k}{h^2}\sin^2(\frac{\beta_y h}{2})}\\
\quad\\
|\xi| \le 1 \iff |1+\frac{2hk\gamma}{h} + \frac{2k}{h^2} - \frac{2k}{h^2} e^{-i\beta_xh}+\frac{4k}{h^2}\sin^2(\frac{\beta_y h}{2})| \ge 1 \\ 
\quad
\\
|1+\frac{2k\gamma}{h} + \frac{2k}{h^2}(1-\cos(\beta_xh)) +\frac{4k}{h^2}\sin^2(\frac{\beta_y h}{2}) +\frac{2ki}{h^2}\sin(\beta_xh)| \ge 1  
\quad\\ 
(1+\frac{2k\gamma}{h} + \frac{2k}{h^2}(1-\cos(\beta_xh)) +\frac{4k}{h^2}\sin^2(\frac{\beta_y h}{2}))^2 +\frac{4k^2}{h^4}\sin^2(\beta_xh) \\
\quad \\
\ge (1+\frac{2k\gamma}{h} )^2 >  1 \Rightarrow |\xi| \le1 \quad \forall \beta_x,\beta_y
$$ 

The same argumentation can be used for the left boundary as well, with the only difference being the sign in the exponential which does not affect the stability. Thus we can conclude our boundary equations gives us Von Neumann stability. 

Lastly we have to show von Neumann stability for the corner equations. The only difference that occurs in the analysis of the different corners is the sign of the exponentials, so we will take into account all the cases in the same computation. In addition to avoid redoing the same computations, we will skip the step of solving for $\xi$ and start off with the following inequality we have shown before. Inserting $U_p^n = \xi^n e^{i\beta0}$ into one of the corner equation yields the following:  
$$
|1 + \frac{4k\gamma}{h} + \frac{k}{h^2}(4 -2e^{\pm i \beta_xh} -2e^{\pm i\beta_yh})| \ge 1 \\ 
\quad \\
|1 + \frac{4k\gamma}{h} + \frac{k}{h^2}(4 -2\cos(\beta_xh)-2\cos(\beta_yh)) +\frac{ik}{h^2}(\mp \sin(\beta_xh) \mp \sin(\beta_yh) )| \\ 
\quad \\ 
(1 + \frac{4k\gamma}{h} + \frac{k}{h^2}(4 -2\cos(\beta_xh)-2\cos(\beta_yh)))^2 +(\frac{k}{h^2}(\mp \sin(\beta_xh) \mp \sin(\beta_yh)))^2 \\  
\quad\\
\ge (1 + \frac{4k\gamma}{h} )^2 > 1 \Rightarrow |\xi| \le1 \quad \forall \beta_x, \beta_y. 
$$

Due to all our equations for for interior, boundary and corners satisfies the Von Neumann stability we can conclude that our scheme is Von Neumann stable.