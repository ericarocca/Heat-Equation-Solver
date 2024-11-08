import numpy as np
import pytest
import configparser
import hypothesis
from hypothesis import strategies as st
from hypothesis import settings
from hypothesis import given
from function import function_temperature, heat_equation_CN, heat_equation_analytical, check_stability, create_matrices, apply_boundary_conditions, stability

config = configparser.ConfigParser()
config.read('configuration.txt')

length = float(config.get('settings', 'length'))
nx_values = list(map(int, config.get('settings', 'nx_values').split(',')))
time = float(config.get('settings', 'time'))
nt_values = list(map(int, config.get('settings', 'nt_values').split(',')))
alpha = float(config.get('settings', 'alpha'))

#test cases from the configuration file.
test_cases = [
    {"length": length, "nx": nx, "time": time, "nt": nt, "alpha": alpha}
    for nx in nx_values for nt in nt_values
]

def test_initial_conditions():
    """
    Test that the initial conditions of the temperature function are respected.
    
    GIVEN: A spatial grid generated by linspace with values ranging from 0 to 1.
    WHEN: Computing the initial temperature using function_temperature.
    THEN: The resulting values should closely match the expected sine function values.

    """
    x = np.linspace(0, 1.0, 20)
    expected_initial =  np.sin(np.pi * x / 1.0)
    actual_initial = np.array([function_temperature(xi, 1.0) for xi in x])
    np.testing.assert_allclose(actual_initial, expected_initial, rtol=1e-3)

@pytest.mark.parametrize("parameters", test_cases)
def test_check_stability(parameters):
    """
    Test that checks that the parameters meet the stability condition.
    If not, skip the specific case.
    """
    length = parameters["length"]
    nx = parameters["nx"]
    time = parameters["time"]
    nt = parameters["nt"]
    alpha = parameters["alpha"]

    stable_combinations = check_stability(length, time, alpha, [nx], [nt])
    
    if len(stable_combinations) == 0:
        pytest.skip(f"Parameters {parameters} do not meet the stability condition.")
    
    assert len(stable_combinations) > 0, f"Parameters {parameters} do not meet the stability condition."

@pytest.mark.parametrize("parameters", test_cases)
def test_boundary_conditions(parameters):
    """
    Test that checks the enforcing of boundary conditions on the w array.
    
    GIVEN a set of parameters for length, spatial steps (nx), time, time steps (nt), and diffusivity coefficient (alpha)
          that meet the stability requirements for the Crank-Nicolson heat equation solver,
    WHEN the heat equation is solved using these parameters,
    THEN the boundary conditions should be enforced on the w temperature array, meaning:
         - the temperature at the first position (index 0) along the rod remains 0 for all time steps.
         - the temperature at the last position (index -1) along the rod remains 0 for all time steps.

    """
    length = parameters["length"]
    nx = parameters["nx"]
    time = parameters["time"]
    nt = parameters["nt"]
    alpha = parameters["alpha"]
    
    try:
        x, w = heat_equation_CN(length, nx, time, nt, alpha, function_temperature)
    except ValueError as e:
        pytest.skip(str(e))
    
    assert np.all(w[0, :] == 0)
    assert np.all(w[-1, :] == 0)

@pytest.mark.parametrize("parameters", test_cases)
def test_matrices_shape(parameters):
    """
    Test that checks that the w matrix has the expected dimension [nx, nt] for the numerical solution.
    
    GIVEN: A set of parameters for the heat equation (length, nx, time, nt, alpha).
    WHEN: Generating the matrices using heat_equation_CN.
    THEN: The output matrix should have the shape (nx, nt).

    """
    length = parameters["length"]
    nx = parameters["nx"]
    time = parameters["time"]
    nt = parameters["nt"]
    alpha = parameters["alpha"]
    
    try:
        x, w = heat_equation_CN(length, nx, time, nt, alpha, function_temperature)
    except ValueError as e:
        pytest.skip(str(e))
        
    assert w.shape == (nx, nt)

@pytest.mark.parametrize("length, time, alpha, nx, nt, expected", [
    (1.0, 1.0, 0.5, 50, 50, False),
    (2.0, 1.0, 0.05, 20, 10, False),
    (1.0, 1.0, 0.25, 100, 400, False),
    (1.0, 1.0, 0.1, 50, 500, True),
    (2.0, 2.0, 0.02, 100, 1000, True), 
    (2.0, 1.5, 0.01, 150, 1000, True)    
])
def test_stability(length, time, alpha, nx, nt, expected):
    """
    Test the stability function for various scenarios.

    GIVEN: specific parameters (length, time, alpha, nx, nt) representing a setup for the heat equation.
    WHEN: checking stability using the stability() function.
    THEN: the result should match the expected stability outcome.

    Parameters:
    - length, time, alpha, nx, nt : float, float, float, int, int
                                  Inputs for the stability function.
    - expected: bool
               Expected boolean outcome indicating stability.
    """
    result = stability(length, time, alpha, nx, nt)
    assert result == expected

def test_create_matrices():
    """
    Test that the create_matrices function returns matrices A and B with correct properties.
    
    GIVEN: A grid size nx and stability coefficient r.
    WHEN: Creating matrices using create_matrices.
    THEN: The matrices A and B should be square, symmetric, and of size (nx, nx).
    """
    nx = 10
    r = 0.25
    A, B = create_matrices(nx, r)

    #check dimensions
    assert A.shape == (nx, nx), "Matrix A has incorrect dimensions"
    assert B.shape == (nx, nx), "Matrix B has incorrect dimensions"
    #check symmetry
    assert np.allclose(A, A.T), "Matrix A is not symmetric"
    assert np.allclose(B, B.T), "Matrix B is not symmetric"

@pytest.mark.parametrize("nx", [5, 10, 15])
def test_apply_boundary_conditions(nx):
    """
    Test that apply_boundary_conditions correctly applies Dirichlet boundary conditions.
    
    GIVEN: A random matrix of size (nx, nx).
    WHEN: Applying Dirichlet boundary conditions using apply_boundary_conditions.
    THEN: The first and last rows should match identity boundary conditions, while other rows remain unaltered.
    """
    matrix = np.random.rand(nx, nx)
    modified_matrix = apply_boundary_conditions(matrix)

    #top row boundary condition
    assert np.array_equal(modified_matrix[0, :], np.eye(nx)[0, :]), "Boundary conditions not applied correctly at the top row"
    #bottom row boundary condition
    assert np.array_equal(modified_matrix[-1, :], np.eye(nx)[-1, :]), "Boundary conditions not applied correctly at the bottom row"
    #other rows remain unaltered
    assert np.array_equal(modified_matrix[1:-1, 1:-1], matrix[1:-1, 1:-1]), "Internal matrix rows modified incorrectly"


@settings(deadline=None)
@given(
    length=st.floats(min_value=1.0, max_value=10.0),  
    nx=st.integers(min_value=25, max_value=250),
    time=st.floats(min_value=1.0, max_value=10.0),
    nt=st.integers(min_value=25, max_value=250),
    alpha=st.floats(min_value=0.01, max_value=0.5)
)
def test_simulation_convergence(length, nx, time, nt, alpha):
    """
    Test that checks both the numerical and analytical solutions.
    First, the stability condition is verified, and then that the
    numerical solution matrix and the analytical one have the same dimensions.
    In the end, it is verified that the numerical solutions match the analytical ones.
    
    GIVEN: Randomized parameters for the heat equation.
    WHEN: Solving the equation numerically and analytically.
    THEN: The solutions should match in shape and values should be close within tolerance.
    """
    deltax = length / (nx - 1)
    deltat = time / (nt - 1)
    r = alpha * deltat / deltax**2

    assume(r < 0.5)
    
    x_num, w_num = heat_equation_CN(length, nx, time, nt, alpha, function_temperature)
    x_ana, w_ana = heat_equation_analytical(length, nx, time, nt, alpha)

    assert w_num.shape == w_ana.shape, f"Shapes of numerical {w_num.shape} and analytical {w_ana.shape} solutions do not match"

    try:
        np.testing.assert_allclose(w_num, w_ana, rtol=1e-3, atol=5e-4)
    except AssertionError as e:
        print(f"Failed case:\nlength={length}, nx={nx}, time={time}, nt={nt}, alpha={alpha}")
        print(f"deltax={deltax}, deltat={deltat}, r={r}")
        raise e
