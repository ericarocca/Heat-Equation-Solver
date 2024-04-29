import numpy as np

#1D heat equation in a rod
#define Crank-Nicolson function

def heatCN(length, nx, time, nt, alpha):
    #nx = spatial steps
    #nt = time steps
    #alpha = diffusivity coefficient of the medium
    
    
    #create space array
    x_grid = np.linspace(0, length, num=nx) 

    #function of temperature
    function_temperature = np.sin( np.pi * x_grid / length)
    #function_temperature = np.cos(np.pi * x_grid)
    
    #create heat array 
    w = np.zeros([nx, nt])
    
    #create matrices for calculation
    A = np.eye(nx)
    B = np.eye(nx)

    #create vector b of boundary conditions
    b = np.zeros(nx)
    
    for i in range(nx):
        w[i,:] = function_temperature[i]
    
    #Dirichlet boundary conditions
    w[0,:] = w[-1,:] = 0
    
    #calculate coefficient
    deltax = length / (nx - 1)
    
    if nt > 1:
        deltat = time / (nt - 1)
        r = alpha * deltat / (deltax) ** 2
        
        #create matrices
        A = np.eye(nx) - ((np.eye(nx, k=1) + np.eye(nx, k=-1) - (1 / r + 2) * np.eye(nx)) * r)
        B = np.eye(nx) + ((np.eye(nx, k=1) + np.eye(nx, k=-1) + (1 / r - 2) * np.eye(nx)) * r)
        
        for i in range(nt):
            #vector b of boundary conditions in RHS
            b[0] = r * 2 * w[0,i]
            b[-1] = r * 2 * w[-1,i]
            
            #initialize d with the same shape as w[:, i]
            d = np.zeros_like(w[:, i])  
            
            #compute matrix-vector product between B and the i-th column of w
            for j in range(nx):
                d[j] = np.dot(B[j], w[:, i])
                
            #add the boundary conditions
            d[0] += b[0]
            d[-1] += b[-1]
            
            
            #solve
            w[:,i] = np.linalg.solve(A, d)
            #impose boundary conditions
            w[0,:] = w[-1,:] = 0
        
    
    #return solution    
    return x_grid, w, b, A, B
