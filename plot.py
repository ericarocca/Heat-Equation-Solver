from function import heat_equation_CN, function_temperature
import numpy as np
import matplotlib.pyplot as plt

def plot_solutions(x, w, wa, nt, time):
    """
    This function plots the comparison between the numerical and the analytical 
    solution of the heat equation evolving through different time steps

    Parameters
    ----------
    x : spatial steps of both the solutions.
    w : array of the temperature calculated with the numerical method.
    wa : array of the temperature calculated with the analytical method.
    nt : time steps of both the solutions.
    time : total time of the evolution.
    """
    plt.figure(figsize=(12, 6))
    timesteps = [0, int(nt/3), int(2*nt/3), nt-1]
    
    for i in timesteps:
        plt.plot(x, w[:, i], label=f'Numerical t={i*time/(nt-1):.2f}')
        plt.plot(x, wa[:, i], '--', label=f'Analytical t={i*time/(nt-1):.2f}')
    
    plt.xlabel('Position')
    plt.ylabel('Temperature')
    plt.legend()
    plt.title('Comparison of Numerical and Analytical Solutions of the Heat Equation')
    plt.show()

def plot_surface_solution(x, w, nt, time):
    """
    This function plots the temperature as a function of both position and time.
    
    Parameters
    ----------
    x : spatial steps of the numerical solution.
    w : array of the temperature calculated with the numerical method.
    nt : time steps of the numerical solution.
    time : total time of the evolution.
    """
    X, T = np.meshgrid(x, np.linspace(0, time, nt))
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, T, w.T, cmap='viridis')
    ax.set_xlabel('Position')
    ax.set_ylabel('Time')
    ax.set_zlabel('Temperature')
    ax.set_title('3D Surface Plot of the Heat Equation Solution')
    plt.show()
