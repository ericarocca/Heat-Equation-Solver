<h1>Heat equation numerical solutions</h1>
The heat equation is a parabolic partial differential equation (also known as PDE), which in one dimension is expressed as:

$\frac{\partial u}{\partial t} = \alpha \frac{\partial^2 u}{\partial x^2}$

where $x$ and $t$ are the spatial and time coordinates, $u$ is the temperature at the point of coordinates $x$ and $\alpha$ is the coefficient of thermal diffusivity of the medium.

The program aims to solve it using the Crank-Nicolson method, an implicit finite difference method.
The discretization of the heat equation using the Crank-Nicolson method is then:

$-r u_{i+1}^{n+1} + (1 + 2r)u_i^{n+1} - r u_{i-1}^{n+1} = r u_{i+1}^n + (1 - 2r) u_i^n + r u_{i-1}^n$

with $r = \frac{\alpha \Delta t}{2(\Delta x)^2}$

<h2>Structure of the project</h2>
