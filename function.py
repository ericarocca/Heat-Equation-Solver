import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def check_stability(length, time, alpha, nx_values, nt_values):
    """
    Checking for the stable parameters.

    Parameters
    ----------
    length : length of the rod.
    time : time of the evolution.
    alpha : diffusivity cooefficient of the medium.
    nx_values : spatial steps.
    nt_values : time steps.

    Returns
    -------
    stable_combinations : array with the combination of parameters that respect the condition of r < 0.5.

    """
    
    stable_combinations = []
    
    for nx in nx_values:
        for nt in nt_values:
            deltax = length / (nx - 1)
            deltat = time / (nt - 1)
            r = alpha * deltat / deltax**2
            if r < 0.5:
                stable_combinations.append((length, time, nx, nt, r))
    
    return stable_combinations
def function_temperature(x, length):
    """
    Generate the initial temperature distribution.
       
    Parameters
        x : spatial steps on the rod.
        length : length of the rod.
    
    Returns:
        The initial temperature distribution.
    """
        
    #return np.sin(np.pi * x)
    return 6 * np.sin(np.pi * x / length)

def heat_equation_CN(length, nx, time, nt, alpha, function_temperature):
    """
    The function calculates the numerical solution of the heat equation using Crank-Nicolson method.
    
    Parameters
    ----------
        length : length of the rod.
        nx : spatial steps.
        time : time of the evolution.
        nt : time steps.
        alpha : diffusivity coefficient of the medium.
        
    Returns
    -------
        x_grid : array of spatial coordinates along the rod with nx points.
        w : array of the temperature calculated with the Crank-Nicolson method of dimensions [nx, nt].
    
    Raises
    ------
        ValueError: If the stability condition is not respected
    """
    
    x_grid = np.linspace(0, length, num=nx)
    w = np.zeros([nx, nt])

    for i in range(nx):
        w[i, 0] = function_temperature(x_grid[i], length)

    w[0, :] = w[-1, :] = 0

    deltax = length / (nx - 1)
    deltat = time / (nt - 1)
    r = alpha * deltat / deltax**2

    #print(f"deltax: {deltax}, deltat: {deltat}, r: {r}")
    
    if r > 0.5:
        raise ValueError(f"The scheme is unstable: r = {r}")

    A = np.eye(nx) - r/2 * (np.eye(nx, k=1) + np.eye(nx, k=-1) - 2 * np.eye(nx))
    B = np.eye(nx) + r/2 * (np.eye(nx, k=1) + np.eye(nx, k=-1) - 2 * np.eye(nx))
    
    #boundary conditions on matrices
    A[0, :] = A[-1, :] = 0
    A[0, 0] = A[-1, -1] = 1
    B[0, :] = B[-1, :] = 0
    B[0, 0] = B[-1, -1] = 1

    for i in range(1, nt):
        d = B @ w[:, i-1]
        d[0] = d[-1] = 0
        w[:, i] = np.linalg.solve(A, d)
        if i == 1 or i == nt - 1:  # Print for the first and last time steps
            print(f"Numerical solution at t={i * deltat}:")
            print(w[:, i])

    
    return x_grid, w

def heat_equation_analytical(length, nx, time, nt, alpha):
    """
    The function calculates the analytical solution of the 1D heat equation.
    
    Parameters
    ----------
    length : length of the rod.
    nx : spatial steps.
    time : evolution time.
    nt : time steps.
    alpha : diffusivity coefficient of the medium.

    Returns
    -------
    x : array of spatial coordinates along the rod with nx points.
    wa : array of the temperature calculated with a Fourier sine series.
    """
    
    wa = np.zeros([nx, nt])
    t = np.linspace(0, time, num=nt)
    x = np.linspace(0, length, num=nx)
    
    for i in range(nt):
        wa[:, i] = 6 * np.sin(np.pi * x / length) * np.exp(-alpha * (np.pi / length)**2 * t[i])
        #wa[:, i] = np.sin(np.pi * x ) * np.exp(-alpha * (np.pi / length)**2 * t[i])
        wa[0, i] = wa[-1, i] = 0
        if i == 1 or i == nt - 1:  # Print for the first and last time steps
            print(f"Analytical solution at t={t[i]}:")
            print(wa[:, i])
    
    return x, wa
