from function import heatCN
import numpy as np
#import pytest

test_cases = [
    #peat thermal diffusivity
    {"length": 1.0, "nx": 2, "time": 1.0, "nt": 2, "alpha": 0.1},
    {"length": 21.0, "nx": 300, "time": 30.0, "nt": 1200, "alpha": 0.1},
    #granite thermal diffusivity
    {"length": 47.0, "nx": 129, "time": 1.0, "nt": 50, "alpha": 1.3}
    ]

def test_boundary_conditions(length, nx, time, nt, alpha):
    
    for params in test_cases:
        length, nx, time, nt, alpha = params
    x_grid, w = heatCN(length, nx, time, nt, alpha)
    
    assert np.all(w[0,:] == 0)
    assert np.all(w[-1,:] == 0)
    
    
def test_matrices_construction():
    for params in test_cases:
        length, nx, time, nt, alpha = params
    x_grid, w, b, _, _, A, B = heatCN(length, nx, time, nt, alpha)
    
    #check dimensions of matrices A, B, C, w
    assert A.shape == (nx, nx)
    assert B.shape == (nx, nx)
    assert w.shape == (nx, nt)
    
    #check dimension of array b
    assert b.shape == (nx,)
        