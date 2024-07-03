import numpy as np
import pytest
import hypothesis
from hypothesis import strategies as st
from hypothesis import settings
from hypothesis import given
from function import function_temperature, heat_equation_CN, heat_equation_analytical

test_cases = [
    # peat thermal diffusivity
    {"length": 1.0, "nx": 2, "time": 1.0, "nt": 2, "alpha": 0.1},
    {"length": 21.0, "nx": 300, "time": 30.0, "nt": 1200, "alpha": 0.1},
    # granite thermal diffusivity
    {"length": 47.0, "nx": 129, "time": 1.0, "nt": 50, "alpha": 1.3}
]

TOLERANCE = 1e-3

def test_initial_conditions():
    x = np.linspace(0, 1.0, 20)
    expected_initial =  np.sin(np.pi * x / 1.0)
    actual_initial = np.array([function_temperature(xi, 1.0) for xi in x])
    np.testing.assert_allclose(actual_initial, expected_initial, rtol=TOLERANCE)

@pytest.mark.parametrize("parameters", test_cases)
def test_boundary_conditions(parameters):
    length = parameters["length"]
    nx = parameters["nx"]
    time = parameters["time"]
    nt = parameters["nt"]
    alpha = parameters["alpha"]

    x_grid, w = heat_equation_CN(length, nx, time, nt, alpha, function_temperature)
    
    assert np.all(w[0, :] == 0)
    assert np.all(w[-1, :] == 0)

@pytest.mark.parametrize("parameters", test_cases)
def test_matrices_shape(parameters):
    length = parameters["length"]
    nx = parameters["nx"]
    time = parameters["time"]
    nt = parameters["nt"]
    alpha = parameters["alpha"]
    
    x_grid, w = heat_equation_CN(length, nx, time, nt, alpha, function_temperature)
        
    assert w.shape == (nx, nt)
