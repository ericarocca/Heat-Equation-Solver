import numpy as np
import math

#1D heat equation in a rod
#define all the functions
#Function of temperature
#Crank-Nicolson function
#Analytical function

def function_temperature(x, length):
    return np.sin(np.pi * x)

def heat_equation_CN(length, nx, time, nt, alpha, function_temperature):
    #nx = spatial steps
    #nt = time steps
    #alpha = diffusivity coefficient of the medium
    
    
    #create space array
    x_grid = np.linspace(0, length, num=nx) 

    #function of temperature
#    function_temperature = np.sin( np.pi * x_grid / length)
#    function_temperature = np.cos(np.pi * x_grid)
#    np.ones_like(x) * 20
    
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




#rivedi k

def heat_equation_analytical(length, nx, time, nt, k):
    
    #create heat array
    wa = np.zeros([nx, nt])
    
    #create time and space arrays
    t = np.linspace(0, time, num=nt)
    x = np.linspace(0, length, num=nx)
    
    #general solution
    #first create the fragments of Fourier series
    #then they combine to form the solution
    for n in range (2,500):
        for i in range (nt-1):
            e = math.exp(-np.pi**2 * k * n**2 *t[i])
            s = np.sin(n * np.pi * x[:])
            #c = (2 * n * ((-1)**n + 1) ) / (np.pi * (n + 1) * (n - 1))
            #s = np.sin((n * np.pi * x[:]))
            #e = math.exp(-0.5*((n*np.pi)**2)*t[i])
            
          
            wa[:,i] = wa[:,i] + (s * e)
    
    
    return x, wa

def analytical_solution(x, timesteps, alpha, L, N):
    x = np.linspace(0, L, N)
    solutions = np.zeros((len(timesteps), len(x)))
    for i, t in enumerate(timesteps):
        solution_t = np.zeros_like(x)
        for n in range(1, 10):  # Adjust range for more terms in the series
            solution_t += np.exp(-np.pi**2 * alpha * n**2 * t) * np.sin(n * np.pi * x / L)
        solutions[i, :] = solution_t
    return solutions


#np.exp(-np.pi**2 * alpha * n**2 * t) * np.sin(n * np.pi * x / L)
