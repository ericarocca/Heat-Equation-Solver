import numpy as np

def stability(length, time, alpha, nx, nt):
    """
    Check if the given parameters are stable for the Crank-Nicolson method.

    Parameters
    ----------
    length : float
        Length of the rod.
    time : float
        Time of the evolution.
    alpha : float
        Diffusivity coefficient of the medium.
    nx : int
        Number of spatial steps.
    nt : int
        Number of time steps.

    Returns
    -------
    bool
        True if the stability condition (r < 0.5) is met, False otherwise.
    """
    r = calculate_r(length, time, nx, nt, alpha)    
    return r < 0.5  #return True if stable

def check_stability(length, time, alpha, nx_values, nt_values):
    """
    Checking for the stable parameters.

    Parameters
    ----------
    length : floar
            length of the rod.
    time : float
          time of the evolution.
    alpha : float
           diffusivity cooefficient of the medium.
    nx_values : list of int
               spatial steps.
    nt_values : list of int
               time steps.

    Returns
    -------
    stable_combinations : array with the combination of parameters that respect the condition of r < 0.5.

    """
    
    stable_combinations = []
    
    for nx in nx_values:
        for nt in nt_values:
            r = calculate_r(length, time, alpha, nx, nt)
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
    return np.sin(np.pi * x / length)

def calculate_r(length, time, nx, nt, alpha):
    """
    Calculate the stability factor r.

    Parameters
    ----------
    length : float
            Length of the rod.
    time : float
          Time of the evolution.
    nx : int
        Number of spatial steps.
    nt : int
        Number of time steps.
    alpha : float
           Diffusivity coefficient of the medium.

    Returns
    -------
    r : float
        Stability factor (alpha * deltat / deltax**2).
    """
    deltax = length / (nx - 1)
    deltat = time / (nt - 1)
    return alpha * deltat / deltax**2


def create_matrices(nx, r):
    """
    Create matrices A and B for the Crank-Nicolson method.

    Parameters
    ----------
    nx : int
        Number of spatial steps.
    r : float
        Stability factor (alpha * deltat / deltax**2).

    Returns
    -------
    A : array
        Matrix A for the Crank-Nicolson method.
    B : array
        Matrix B for the Crank-Nicolson method.
    """
    A = np.eye(nx) - r/2 * (np.eye(nx, k=1) + np.eye(nx, k=-1) - 2 * np.eye(nx))
    B = np.eye(nx) + r/2 * (np.eye(nx, k=1) + np.eye(nx, k=-1) - 2 * np.eye(nx))
    
    return A, B

def apply_boundary_conditions(matrix):
    """
    Apply Dirichlet boundary conditions to the matrix for the Crank-Nicolson method.

    Parameters
    ----------
    matrix : array
            Matrix for the Crank-Nicolson method.

    Returns
    -------
    matrix : array
            Modified matrix with Dirichlet boundary conditions applied.
    """
    matrix[0, :] = matrix[-1, :] = 0
    matrix[0, 0] = matrix[-1, -1] = 1

    return matrix

def heat_equation_CN(length, nx, time, nt, alpha, function_temperature):
    """
    The function calculates the numerical solution of the heat equation using Crank-Nicolson method.
    
    Parameters
    ----------
        length : float
                length of the rod.
        nx : int
            spatial steps.
        time : float
              evolution time.
        nt : int
            time steps.
        alpha : float
               diffusivity coefficient of the medium.
        
    Returns
    -------
        x : array
           spatial coordinates along the rod with nx points.
        w : array
           temperature calculated with the Crank-Nicolson method, dimensions [nx, nt].
    
    Raises
    ------
        ValueError: If the stability condition is not respected
    """

    if not stability(length, time, alpha, nx, nt):
        raise ValueError("The scheme is unstable for the provided parameters.")

    x = np.linspace(0, length, num=nx)
    w = np.zeros([nx, nt])

    for i in range(nx):
        w[i, 0] = function_temperature(x[i], length)

    w[0, :] = w[-1, :] = 0

    r = calculate_r(length, time, nx, nt, alpha)
    
    A, B = create_matrices(nx, r)
    
    A = apply_boundary_conditions(A)
    B = apply_boundary_conditions(B)

    for i in range(1, nt):
        d = B @ w[:, i-1]
        d[0] = d[-1] = 0
        w[:, i] = np.linalg.solve(A, d)

    return x, w

def heat_equation_analytical(length, nx, time, nt, alpha):
    """
    The function calculates the analytical solution of the 1D heat equation.
    
    Parameters
    ----------
        length : float
                length of the rod.
        nx : int
            spatial steps.
        time : float
             evolution time.
        nt : int
            time steps.
        alpha : float
               diffusivity coefficient of the medium.

    Returns
    -------
    x : array
       spatial coordinates along the rod with nx points.
    wa : array
        temperature calculated with the Fourier sine series.
    """
    
    stable_combinations = check_stability(length, time, alpha, [nx], [nt])
    if not stable_combinations:
        raise ValueError("The scheme is unstable for the provided parameters.")

    wa = np.zeros([nx, nt])
    t = np.linspace(0, time, num=nt)
    x = np.linspace(0, length, num=nx)
    
    for i in range(nt):
        wa[:, i] = np.sin(np.pi * x / length) * np.exp(-alpha * (np.pi / length)**2 * t[i])
        wa[0, i] = wa[-1, i] = 0
        
    return x, wa
