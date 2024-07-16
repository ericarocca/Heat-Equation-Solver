<h1>Heat equation numerical solutions</h1>
The heat equation is a parabolic partial differential equation (PDE) that shows how heat changes over time in a solid medium.
In one dimension is expressed as:

$\frac{\partial w}{\partial t} = \alpha \frac{\partial^2 w}{\partial x^2}$

where $x$ and $t$ are the spatial and time coordinates, $u$ is the temperature at the point of coordinates $x$ and $\alpha$ is the coefficient of thermal diffusivity of the medium.

The program aims to examine the heat diffusion in a one dimension rod with fixed boundary conditions at both ends and an initial temperature distribution given by a sinusoidal function. The equation is solved numerically using the Crank-Nicolson method, an implicit finite difference method, and it is then compared with the analytical solution obtained with a Fourier sine series.\
The discretization of the heat equation using the Crank-Nicolson method is then:

$-r w_{i+1}^{n+1} + (1 + 2r)w_i^{n+1} - r w_{i-1}^{n+1} = r w_{i+1}^n + (1 - 2r) w_i^n + r w_{i-1}^n$

with $r = \frac{\alpha \Delta t}{2(\Delta x)^2}$.\
The equation can also be written in matrix form as:

$w_{n+1} = A^{-1}Bw_n$

where $A$ and $B$ are tridiagonal matrices

<h2>Structure of the project</h2>

There are four blocks in this project:
* In the [function.py](https://github.com/ericarocca/Heat-Equation-Solver/blob/main/function.py) file there are all the functions.
* In the [test.py](https://github.com/ericarocca/Heat-Equation-Solver/blob/main/test.py) file all the functions are tested so that all of them work properly as well as the program.
* In the [plot.py](https://github.com/ericarocca/Heat-Equation-Solver/blob/main/plot.py) file there are two functions, one that plot the comparison between the numerical and the analytical solution of the heat equation and the other that show the surface plot of the numerical solution.
* In the [simulation.py](https://github.com/ericarocca/Heat-Equation-Solver/blob/main/simulation.py) file there is the main part of the code, where the numerical and analytical solutions are calculated and then plotted
