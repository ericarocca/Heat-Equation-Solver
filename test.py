from function import heatCN
import pytest
import numpy as np


test_cases = [
    #peat thermal diffusivity
    {"length": 1.0, "nx": 2, "time": 1.0, "nt": 2, "alpha": 0.1},
    {"length": 21.0, "nx": 300, "time": 30.0, "nt": 1200, "alpha": 0.1},
    #granite thermal diffusivity
    {"length": 47.0, "nx": 129, "time": 1.0, "nt": 50, "alpha": 1.3}
    ]


@pytest.mark.parametrize("parameters", test_cases)
def test_boundary_conditions(parameters):
    length = parameters["length"]
    nx = parameters["nx"]
    time = parameters["time"]
    nt = parameters["nt"]
    alpha = parameters["alpha"]
    
    x_grid, w, b, A, B = heatCN(length, nx, time, nt, alpha)
    
    assert np.all(w[0,:] == 0)
    assert np.all(w[-1,:] == 0)    
    

@pytest.mark.parametrize("parameters", test_cases)
def test_matrices_shape(parameters):
    
    length = parameters["length"]
    nx = parameters["nx"]
    time = parameters["time"]
    nt = parameters["nt"]
    alpha = parameters["alpha"]
   
    length = float(length)
    nx = int(nx)
    time = float(time)
    nt = int(nt)
    alpha = float(alpha)
        
    x_grid, w, b, A, B = heatCN(length, nx, time, nt, alpha)
        
    assert A.shape == (nx, nx)
    assert B.shape == (nx, nx)
    assert w.shape == (nx, nt)
    assert b.shape == (nx,)

    




if __name__ == "main":
    pass        
