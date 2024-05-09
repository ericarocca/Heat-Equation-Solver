from function import heat_equation_CN, function_temperature
import numpy as np
import matplotlib.pyplot as plt

#Convergence study
length = 1.0
time = 1.0
alpha = 2.73  #quartzite typical diffusivity coefficient

#different time step sizes
nt_values = [ 55, 100, 250, 500]

#different spatial step sizes
nx_values = [ 55, 100, 250, 500]


fig, axes = plt.subplots(len(nt_values), len(nx_values), figsize=(12, 12), sharex=True, sharey=True)

for i, nt in enumerate(nt_values):
    for j, nx in enumerate(nx_values):
        x_grid = np.linspace(0, length, num=nx)
        x_grid, w, b, A, B = heat_equation_CN(length, nx, time, nt, alpha, function_temperature(x_grid))
        axes[i, j].plot(x_grid, w[:, -1], label=f'nx={nx}\nnt={nt}')
        axes[i, j].set_title(f'nx={nx}, nt={nt}')
        axes[i, j].legend()

fig.text(0.5, 0.04, 'x', ha='center')
fig.text(0.04, 0.5, 'Temperature', va='center', rotation='vertical')
plt.tight_layout()
plt.show()
