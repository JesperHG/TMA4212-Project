The spacial error term will be dominated mostly by the use of false boundaries to solve the Robin boundary conditions. Rescaling the boundary conditions gave us the equation $u_{n^*} = -\frac{hL}{k}u$. Define $\gamma := \frac{hL}{k}u$, and consider the right bound of the domain, giving us $u_{n^*} = u_{x^*}$. For simplicity's sake we'll omit the use of * to denote the rescaled parameters, but we're still looking at the rescaled system. In an arbitrary right boundary point we then have $u_x = -\gamma u$, and our standard 5-point stencil for computing $u$:

$\Delta_h u = \frac{1}{h^2}(u_N + u_S + u_E + u_W - 4u).$ Since we're at the right boundary, we need to approximate $u_E$ using the Robin condition. We can approximate our boundary equation with a central difference method:
$$\frac{1}{2h}(u_E - u_W) = -\gamma u \\
\Rightarrow u_E = u_W - 2h\gamma u = u_W + 2h u_x.$$ All that remains is the Taylor expansion: 
$u_N = u + hu_y+\frac{h^2}{2}u_{yy} + \frac{h^3}{6}u_{3y} + \mathcal{O}(h^4) \\ 
u_S = u - hu_y+\frac{h^2}{2}u_{yy} - \frac{h^3}{6}u_{3y} + \mathcal{O}(h^4) \\ 
u_W = u - hu_x+\frac{h^2}{2}u_{xx} - \frac{h^3}{6}u_{3x} + \mathcal{O}(h^4).$ We see that $u_N + u_S = 2u + h^2u_{yy} + \mathcal{O}(h^4)$. Inserting all this into our approximation for the laplacian on the boundary gives us $$\Delta_h u = \frac{1}{h^2}(u_N + u_S + u_E + u_W - 4u) \\ = \Delta_h u = \frac{1}{h^2}(u_N + u_S  + 2u_W +2hu_x- 4u)\\ = \frac{1}{h^2}(h^2u_{yy} + 2(u-hu_x+\frac{h^2}{2}u_{xx}-\frac{h^3}{6}u_{3x}) + 2hu_x - 2u + \mathcal{O}(h^4)) \\ = \Delta u - \frac{h}{3}u_{3x} + \mathcal{O}(h^2).$$

So on the boundaries the local truncation error is of degree $1$, while we have degree $2$ on the interior points. 