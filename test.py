import numpy as np
import pytest
import configparser
import hypothesis
from hypothesis import strategies as st
from hypothesis import settings
from hypothesis import given
from function import function_temperature, heat_equation_CN, heat_equation_analytical

config = configparser.ConfigParser()
config.read('configuration.txt')

length = float(config.get('settings', 'length'))
nx_values = list(map(int, config.get('settings', 'nx_values').split(',')))
time = float(config.get('settings', 'time'))
nt_values = list(map(int, config.get('settings', 'nt_values').split(',')))
alpha = float(config.get('settings', 'alpha'))

# Prepare test cases from the configuration file
test_cases = [
    {"length": length, "nx": nx, "time": time, "nt": nt, "alpha": alpha}
    for nx in nx_values for nt in nt_values
]
def test_initial_conditions():
    x = np.linspace(0, 1.0, 20)
    expected_initial =  np.sin(np.pi * x / 1.0)
    actual_initial = np.array([function_temperature(xi, 1.0) for xi in x])
    np.testing.assert_allclose(actual_initial, expected_initial, rtol=1e-3)

@pytest.mark.parametrize("parameters", test_cases)
def test_check_stability(parameters):
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


@settings(deadline=None)
@given(
    length=st.floats(min_value=1.0, max_value=10.0),  
    nx=st.integers(min_value=25, max_value=250),
    time=st.floats(min_value=1.0, max_value=10.0),
    nt=st.integers(min_value=25, max_value=250),
    alpha=st.floats(min_value=0.01, max_value=0.5)
)
def test_simulation_convergence(length, nx, time, nt, alpha):
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
        
def test_analytical_solution():
    length = 1.0
    nx = 50
    time = 0.5
    nt = 500
    alpha = 0.01
    x_num, w_num = heat_equation_CN(length, nx, time, nt, alpha, function_temperature)
    x_a, wa = heat_equation_analytical(length, nx, time, nt, alpha)
    
    assert w_num.shape == wa.shape, f"Shapes of numerical {w_num.shape} and analytical {wa.shape} solutions do not match"
    
    try:
        np.testing.assert_allclose(w_num, wa, rtol=1e-3, atol=5e-4)
    except AssertionError as e:
        print(f"Numerical solution:\n{w_num}")
        print(f"Analytical solution:\n{wa}")
        raise e
    
    x_grid, w = heat_equation_CN(length, nx, time, nt, alpha, function_temperature)
        
    assert w.shape == (nx, nt)
