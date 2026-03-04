---
output:
  pdf_document: default
  html_document: default
---
## Mathematical model 

We have the following model:
$$
\begin{align*}
\rho c_p \frac{\partial T}{\partial t} &= k\Delta T + \rho H_r \frac{\partial \alpha}{\partial t} \quad \text{on} \; \Omega, \\
-k\frac{\partial T}{\partial n} &= h(T-T_m) \quad\quad\quad \text{on} \;\partial\Omega.
\end{align*}
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
\end{align*}v
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

We plug this into the model on $\Omega$, and get:
$$
\begin{align*}
\rho c_p\frac{\partial T}{\partial t} &= k\Delta T + \rho H_r \frac{\partial \alpha}{\partial t} \\
\rho c_p \frac{T_m-T_0}{t_\text{diff}}u_{t^*} &= k\frac{T_m-T_0}{L^2}\Delta^*u + \rho H_r \frac{1}{t_\text{diff}}\alpha_{t^*}\\ 
\rho c_p \frac{T_m-T_0}{t_\text{diff}}u_{t^*} &= \rho c_p\frac{T_m-T_0}{t_\text{diff}}\Delta^*u + \rho H_r \frac{1}{t_\text{diff}}\alpha_{t^*}\\
u_{t^*} &= \Delta^*u + \frac{H_r}{c_p(T_m-T_0)} \alpha_{t^*}\\  
u_{t^*} &= \Delta^*u + c \alpha_{t^*}, \quad c= \frac{H_r}{c_p(T_m-T_0)}. 
\end{align*}
$$

This equation is dimensionless, which is what we wanted. We do the same for the model on $\partial\Omega$, and get:
$$
\begin{align*}
-k\frac{\partial T}{\partial n} &= h(T-T_m) \\
-k\frac{T_m - T_0}{L}u_{n^*} &= h(T-T_m) \\
u_{n^*} &= -\frac{hL}{k}u \\
u_{n^*} + \gamma u &= 0, \quad \gamma=\frac{hL}{k}.
\end{align*}
$$
Now we have the model in dimensionless form. 
 
## Numerical model 
Let the grid be defined as: 
$$
\begin{align*}
G = [x_i,y_j,t_n]\\ 
x_i = ih, \quad y_j = jh, \quad i,j = 1,2,...,M \quad h= \frac{1}{M}\\
t_n = nk, \quad n= 1,2,...,N \ \quad k = \frac{1}{N} 
\end{align*}
$$

We use the following approximations for dicretizations with known error term: 
$$
\begin{align*}
\frac{1}{k}\nabla_k U_p =  \frac{U_p^n - U_p^{n-1}}{k} = u_t -\frac{k}{2}u_{xx}(\xi_m)\\
\frac{1}{h^2}(\delta_x^2 + \delta_y^2) = \Delta_* u +\frac{h^2}{12}(u_{4x}(\xi_m)+u_{4y}(\xi_m))
\end{align*}
$$

Inserting the discretizations in to our dimensionless interior equation yields the following: 
$$
\begin{align*} 
\frac{U_p^{n} - U_p^{n-1}}{k} = \frac{1}{h^2}(-4U_P^n+U_W^n+U_E^n+U_S^n+U_N^n) +kc\alpha_{t^*} \\
U_p^{n} = U_p^{n-1} + \frac{k}{h^2}(-4U_P^n+U_W^n+U_E^n+U_S^n+U_N^n) +kc \alpha_{t^*} \\ 
(4+\frac{h^2}{k})U_p^n -U_W^n -U_N^n-U_S^n-U_E^n = \frac{h^2}{k}U_p^{n-1} + h^2c\alpha_{t^*}\\ 
\end{align*}  
$$



Which gives us the following linear system for the interior points:
$$
\begin{align*}
A_hU_p^n = b\\
A = \text{tridiag}(-I_{M},B,-I_{M}),\quad B = \text{tridiag}(-1,4+\frac{h^2}{k},-1)\\
b = \frac{h^2}{k}U_p^{n-1} + h^2c \alpha_{t^*}\\ 
\end{align*}  
$$

For the boundary equation we use the following approximation: 
$$
\begin{align*}
u_{n^*} + \gamma u = 0 \\ 
\frac{1}{2h}\delta_0u = u_{n^*} + \frac{h^2}{6} u_{3n^*}\\ 
\end{align*}
$$

This yields the following equations for the different boundaries: 
$$
\begin{align*}
\frac{1}{2h}(u_E-u_W) &= -\gamma u \quad \text{Right boundary}\\
\quad\\
\frac{1}{2h}(u_W-u_E) &= -\gamma u \quad \text{Left boundary}\\ 
\quad\\
\frac{1}{2h}(u_N-u_S) &= -\gamma u \quad\text{Upper boundary}\\
\quad\\ 
\frac{1}{2h}(u_S-u_N) &= -\gamma u \quad\text{Bottom boundary}\\
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
-1 & a & -1 & \cdots & 0 \\
0 & -1 & a & \ddots & \vdots \\
\vdots & \vdots & \ddots & \ddots & -1 \\
0 & 0 & \cdots & -2 & d 
\end{bmatrix}
\\
a &= 4 + \frac{h^2}{k}\\
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

 
## Stability 
First we peform stability analysis on the interior equations:  
$$
\begin{align*}
U_p^{n} = U_p^{n-1} + \frac{k}{h^2}(-4U_P^n+U_W^n+U_E^n+U_S^n+U_N^n) \\ 
\xi^n = \xi^{n-1} + \frac{k}{h^2}\xi^n(-4 + e^{i\beta_xh} +  e^{-i\beta_xh} + e^{i\beta_yh} +  e^{-i\beta_yh}) \xi= 1 + \frac{2k}{h^2}(-2\sin^2(\frac{\beta_xh}{2})-2\sin^2(\frac{\beta_yh}{2}))\\ 
\quad \\ \xi = \frac{1}{1 + \frac{4k}{h^2}(2\sin^2(\frac{b_xh}{2})+2\sin^2(\frac{b_yh}{2}))} \\ \quad\\ |\xi| \le 1 \quad  \forall \beta_x,\beta_y. 
\end{align*} 
$$

Next is the boundary equations, we start by considering the right boundary equation:
$$
\begin{align*}
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
\end{align*}  
$$

The same argumentation can be used for the left boundary as well, with the only difference being the sign in the exponential which does not affect the stability. Thus we can conclude our boundary equations gives us Von Neumann stability. 

Lastly we have to show von Neumann stability for the corner equations. The only difference that occurs in the analysis of the different corners is the sign of the exponentials, so we will take into account all the cases in the same computation. In addition to avoid redoing the same computations, we will skip the step of solving for $\xi$ and start off with the following inequality we have shown before. 
$$
\begin{align*}
|1 + \frac{4k\gamma}{h} + \frac{k}{h^2}(4 -2e^{\pm i \beta_xh} -2e^{\pm i\beta_yh})| \ge 1 \\ 
\quad \\
|1 + \frac{4k\gamma}{h} + \frac{k}{h^2}(4 -2\cos(\beta_xh)-2\cos(\beta_yh)) +\frac{ik}{h^2}(\mp \sin(\beta_xh) \mp \sin(\beta_yh) )| \\ 
\quad \\ 
(1 + \frac{4k\gamma}{h} + \frac{k}{h^2}(4 -2\cos(\beta_xh)-2\cos(\beta_yh)))^2 +(\frac{k}{h^2}(\mp \sin(\beta_xh) \mp \sin(\beta_yh)))^2 \\  
\quad\\
\ge (1 + \frac{4k\gamma}{h} )^2 > 1 \Rightarrow |\xi| \le1 \quad \forall \beta_x, \beta_y. 
\end{align*}
$$
Due to all our equations for for interior, boundary and corners satisfies the Von Neumann stability we can conclude that our scheme is Von Neumann stable. 

## Consistency analysis 
Truncation error for the interior points is given by: 
$$
\begin{align*}
\tau = (L-L_h)u_p = (\partial t-\frac{1}{k}\nabla_k +\frac{1}{h^2}(\delta_x^2 + \delta_y^2)-\Delta)u_p\\ 
\tau = \frac{k}{2}u_{tt}(\lambda_n) + \frac{h^2}{12}(u_{4x}(\xi_m)+u_{4y}(\eta_m))\\
\tau \le O(k + h^2)\quad, \lambda_n \in [t_n,t_{n+1}], \xi_m \in[x_m,x_{m+1}], \eta_m \in [y_m,y_{m+1}]
\end{align*}
$$
Thus our interior equation has order 1 in time and order 2 in space. 

The spacial error term will be dominated mostly by the use of false boundaries to solve the Robin boundary conditions. Rescaling the boundary conditions gave us the equation $u_{n^*} = -\frac{hL}{k}u$. Define $\gamma := \frac{hL}{k}u$, and consider the right bound of the domain, giving us $u_{n^*} = u_{x^*}$. For simplicity's sake we'll omit the use of * to denote the rescaled parameters, but we're still looking at the rescaled system. In an arbitrary right boundary point we then have $u_x = -\gamma u$, and our standard 5-point stencil for computing $u$:

$\Delta_h u = \frac{1}{h^2}(u_N + u_S + u_E + u_W - 4u).$ Since we're at the right boundary, we need to approximate $u_E$ using the Robin condition. We can approximate our boundary equation with a central difference method:
$$\frac{1}{2h}(u_E - u_W) = -\gamma u \\
\Rightarrow u_E = u_W - 2h\gamma u = u_W + 2h u_x.$$ All that remains is the Taylor expansion: 
$$
\begin{align*}
u_N = u + hu_y+\frac{h^2}{2}u_{yy} + \frac{h^3}{6}u_{3y} + \mathcal{O}(h^4) \\ 
u_S = u - hu_y+\frac{h^2}{2}u_{yy} - \frac{h^3}{6}u_{3y} + \mathcal{O}(h^4) \\ 
u_W = u - hu_x+\frac{h^2}{2}u_{xx} - \frac{h^3}{6}u_{3x} + \mathcal{O}(h^4). 
\end{align*}  
$$
We see that $u_N + u_S = 2u + h^2u_{yy} + \mathcal{O}(h^4)$. Inserting all this into our approximation for the laplacian on the boundary gives us:  
$$
\begin{align*}
\Delta_h u = \frac{1}{h^2}(u_N + u_S + u_E + u_W - 4u) \\ = \Delta_h u = \frac{1}{h^2}(u_N + u_S  + 2u_W +2hu_x- 4u)\\ = \frac{1}{h^2}(h^2u_{yy} + 2(u-hu_x+\frac{h^2}{2}u_{xx}-\frac{h^3}{6}u_{3x}) + 2hu_x - 2u + \mathcal{O}(h^4)) \\ = \Delta u - \frac{h}{3}u_{3x} + \mathcal{O}(h^2).
\end{align*}
$$  


The local truncation error will behave similarily on all boundaries, the only differences being if we differentiate with respect to $x$ or $y$, and the sign in front of the first order term. On the corners we will see similar errors, f.ex the top right corner: $$\tau = -\frac{h}{2}(u_{3x} + u_{3y}) + \mathcal{O}(h^2).$$

So on the boundaries the local truncation error is of degree $1$, while we have degree $2$ on the interior points.  



Lastly we have to adress how we handle the source term in the heat equation. We use the following splitting scheme: 

First step is solving the homogeneous heat equation by solving the linear equation: 
$$
A_hU_p^n =  \frac{h^2}{k}U_p^{n-1}  
$$  

Then we use the temperature obtained by the first step to update the sourceterm and use it to solve the inhomogeneous heat equation and updating alpha with a forward difference: 
$$
\begin{align*}
AU_p^{n} &= U_p^{n-1} + h^2c \alpha_{t^*}\\
\alpha^n &= \alpha^{n-1} + k \alpha_{t^*} 
\end{align*}
$$ 


## Discretization in 3D
We used practically the same procedure for the 3D case as for the 2D case of our vulcanization solver. We go from a 5-point stencil in 2D to a 7-point stencil in 3D. We get the following discretization of our problem:

$$
\begin{align*} 
(6+\frac{h^2}{k})U_P^n - U_E^n - U_W^n - U_N^n - U_S^n - U_U^n - U_D^n &= \frac {h^2}{k}U_P^{n-1} + h^2c\alpha_{t^*} \quad \text{Interior}\\
\frac{1}{2h}(u_Q - u_{Q^{-1}}) &= -\gamma u  \quad \text{Boundaries}
\end{align*}
$$

Here we define $Q$ as the ghost point we need to replace, and $Q^{-1}$ as the inner point in the other direction. As an example on the boundary $x=0$, we have:

$$
\frac{1}{2h}(u_W - u_E) = -\gamma u,
$$

and on the boundary $x=L$, we have:

$$
\frac{1}{2h}(u_E - u_W) = -\gamma u.
$$

We use these equations to set up the matrix $A_h$ for the same problem as before:
$$
A_hU_p^n = b = \frac{h^2}{k}U_p^{n-1} + h^2c\alpha_{t^*}.
$$

### Level 1: $M \times M$ tridiagonal blocks

**Inner** (Core is inner points, boundaries are side points):

$$
T_{\text{inner}} = \begin{pmatrix}
6+\frac{h^2}{k}+2\gamma h & -2 \\
-1 & 6+\frac{h^2}{k} & -1 \\
& \ddots & \ddots & \ddots \\
& & -1 & 6+\frac{h^2}{k} & -1 \\
& & & -2 & 6+\frac{h^2}{k}+2\gamma h
\end{pmatrix}
$$

**Side** (Core is side points, boundaries are edge points):

$$
T_{\text{side}} = \begin{pmatrix}
6+\frac{h^2}{k}+4\gamma h & -2 \\
-1 & 6+\frac{h^2}{k}+2\gamma h & -1 \\
& \ddots & \ddots & \ddots \\
& & -1 & 6+\frac{h^2}{k}+2\gamma h & -1 \\
& & & -2 & 6+\frac{h^2}{k}+4\gamma h
\end{pmatrix}
$$

**Edge** (Core is edge points, boundaries are corner points):

$$
T_{\text{edge}} = \begin{pmatrix}
6+\frac{h^2}{k}+6\gamma h & -2 \\
-1 & 6+\frac{h^2}{k}+4\gamma h & -1 \\
& \ddots & \ddots & \ddots \\
& & -1 & 6+\frac{h^2}{k}+4\gamma h & -1 \\
& & & -2 & 6+\frac{h^2}{k}+6\gamma h
\end{pmatrix}
$$

$+0\gamma h$ are inner points, $+2\gamma h$ are side points,  $+4\gamma h$ are edge points, $+6\gamma h$ are corner points.

### Level 2: $M^2 \times M^2$ block-tridiagonal matrices

**B1** (Core is inner matrices, boundaries are side matrices):

$$
B_{1} = \begin{pmatrix}
T_{\text{side}} & -2I \\
-I & T_{\text{inner}} & -I \\
& \ddots & \ddots & \ddots \\
& & -I & T_{\text{inner}} & -I \\
& & & -2I & T_{\text{side}}
\end{pmatrix}
$$

**B2** (Core is side matrices, boundaries are edge matrices):

$$
B_{2} = \begin{pmatrix}
T_{\text{edge}} & -2I \\
-I & T_{\text{side}} & -I \\
& \ddots & \ddots & \ddots \\
& & -I & T_{\text{side}} & -I \\
& & & -2I & T_{\text{edge}}
\end{pmatrix}
$$

Here $I$ is the $M \times M$ identity matrix.

### Level 3: The full $M^3 \times M^3$ system matrix
Core is B1 matrices, boundaries are B2 matrices
$$
A_h = \begin{pmatrix}
B_{2} & -2I \\
-I & B_{1} & -I \\
& \ddots & \ddots & \ddots \\
& & -I & B_{1} & -I \\
& & & -2I & B_{2}
\end{pmatrix}
$$

where $I$ now is the $M^2 \times M^2$ identity matrix.



## Theory on linear curing law

With a linear curing law $$\frac{\partial\alpha}{\partial t} = \begin{cases} 0 \hspace{2.7cm} T \le T_s \\ \lambda (T - T_s) - \gamma \alpha \hspace{0.5cm} T > T_s\end{cases}$$ we can write our heat equation into a linear system. We introduce $\theta = T - T_s$, and we rewrite our equation
$$\rho c_p \frac{\partial T}{\partial t} = k\Delta T + \rho H_r \frac{\partial \alpha}{\partial t} \\ \Rightarrow \frac{\partial \theta}{\partial t} = K\Delta\theta + C(\lambda\theta - \gamma \alpha),$$ with $K = \frac{k}{\rho c_p}$ and $C = \frac{H_r}{c_p}$. So for $T > T_s$ we get the linear system $$\frac{\partial}{\partial t}\begin{pmatrix} \theta \\ \alpha \end{pmatrix} = \begin{pmatrix} K\Delta + C\lambda & -C\gamma \\ \lambda & -\gamma \end{pmatrix}\begin{pmatrix} \theta \\ \alpha \end{pmatrix}.$$

With central differences and using the Forward Euler we can find as expected a local truncation error of $\tau = \mathcal{O}(\tau + \delta^2)$ with $\tau$ being the time step and $\delta$ the space step. What's more interesting is looking at a Von Neumann analysis of the system. Setting $$\theta^n_m = \hat{\theta}\xi^ne^{i\beta p_m}\\ \alpha^n_m = \hat{\alpha}\xi^ne^{i\beta p_m}$$ we get $$\hat{\theta} (\xi^{n+1} - \xi^n) = \frac{\tau}{\delta^2}K \hat{\theta} \xi^n(-4\sin^2(\frac{\beta_x \delta}{2}) -4\sin^2(\frac{\beta_y \delta}{2})) + \tau C \xi^n(\lambda\hat{\theta} - \gamma \hat{\alpha}).$$ Here we use what we already computed in the previous Von Neumann analysis. Dividing by $\xi^n$ and combining the sinus terms into a constant we get $$\xi \hat{\theta} = \hat{\theta} - \frac{\tau}{\delta^2} \mu\hat{\theta} + \tau C( \lambda\hat{\theta} - \gamma\hat{\alpha}).$$ Similarly for $\alpha$: $$\hat{\alpha}(\xi^{n+1} - \xi^n) = \tau\xi^n(\lambda\hat{\theta} - \gamma\hat{\alpha}) \\ \Rightarrow \xi \hat{\alpha} = \hat{\alpha} + \tau(\lambda\hat{\theta} - \gamma\hat{\alpha})$$ giving us the linear system $$\xi \begin{pmatrix} \hat{\theta} \\ \hat{\alpha} \end{pmatrix} = \begin{pmatrix} 1 - \frac{\tau}{\delta^2} \mu + \tau C\lambda & -\tau C \gamma \\ \tau \lambda & 1 - \tau\gamma \end{pmatrix} \begin{pmatrix} \hat{\theta} \\ \hat{\alpha} \end{pmatrix}.$$ For this system to be stable the matrix $A$ must fulfill $\rho(A) < 1$ to be stable. Calculating the eigenvalues of this matrix and finding the values which fulfill this condition is a bit tricky, since we have $4$ free variables we can change: $\tau, \delta,\gamma \text{ and } \lambda$. Instead we'll use Gershgorin's theorem to find the disks in which our eigenvalues must lie in. $$ D_1 = \{z \mid \vert 1 - \frac{\tau}{\delta^2} \mu + \tau C\lambda - z \vert \le \vert\tau C \gamma\vert\} \\ D_2 = \{z \mid \vert 1 -\tau\gamma - z \vert \le \vert\tau\lambda\vert\}.$$ 
Since we want small step sizes for smaller truncation error, the values of $\gamma$ and $\lambda$ will affect the domain of the disks more than the step sizes. Here we see that for $D_1$ to be contain few values larger than $1$, we need $\lambda$ to be negative and have a large absolute value, and to minimize the radius of the disk, $\gamma$ should be small. For $D_2$ we see the opposite, to be centered away from $1$ we need large $\gamma$, and a small radius leads to small $\lambda$. From this we can see that the stability region of the linear system is quite restrictive. Do note that this is just an approximation of the eigenvalues, since we are using Gershgorin's theorem. Even though the disks contain values with absolute value larger than $1$, this doesn't necessarily mean that the eigenvalues $\vert\sigma\vert \ge 1$.



