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
To start the program the user need to follow these steps:
 

There are five blocks in this project:
* In the [configuration.txt](https://github.com/ericarocca/Heat-Equation-Solver/blob/main/configuration.txt) file there are all the parameters used in the [simulation.py](https://github.com/ericarocca/Heat-Equation-Solver/blob/main/simulation.py). For both nx_values and nt_values there is a list of different parameters so that it is possible to verify more than one combination per execution. There are also local paths for saving the solutions array.
* In the [function.py](https://github.com/ericarocca/Heat-Equation-Solver/blob/main/function.py) file there are all the functions used to calculate the solutions of the heat equation.
* In the [test.py](https://github.com/ericarocca/Heat-Equation-Solver/blob/main/test.py) file all the functions are tested so that all of them work properly as well as the program itself.
* In the [plot.py](https://github.com/ericarocca/Heat-Equation-Solver/blob/main/plot.py) file there are two functions, one that plot the comparison between the numerical and the analytical solution of the heat equation and the other that show the surface plot of the numerical solution through time.
* In the [simulation.py](https://github.com/ericarocca/Heat-Equation-Solver/blob/main/simulation.py) file there is the main part of the code, where the numerical and analytical solutions  matrices are calculated, saved on the appropriate path and then plotted.

Some examples of the obtainable results:
![Plot](./Plot/Figure1.png)
![Plot](./Plot/Figure2.png)
![Plot](./Plot/Figure3.png)
![Plot](./Plot/Figure4.png)
