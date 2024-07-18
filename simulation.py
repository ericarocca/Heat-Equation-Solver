import configparser
from function import heat_equation_CN, heat_equation_analytical, function_temperature, plot_solutions, plot_surface_solution, check_stability

config = configparser.ConfigParser()
config.read('configuration.txt')

length = float(config.get('settings', 'length'))
nx_values = list(map(int, config.get('settings', 'nx_values').split(',')))
time = float(config.get('settings', 'time'))
nt_values = list(map(int, config.get('settings', 'nt_values').split(',')))
alpha = float(config.get('settings', 'alpha'))

#Verify the presence of stable combinations and then solve and plot for those.
stable_combinations = check_stability(length, time, alpha, nx_values, nt_values)
if not stable_combinations:
    raise ValueError("No stable combinations found for the provided parameters.")

for combination in stable_combinations:
    chosen_length, chosen_time, chosen_nx, chosen_nt, chosen_r = combination

    print(f"Running simulation with nx={chosen_nx}, nt={chosen_nt}, r={chosen_r}")

    x, w = heat_equation_CN(chosen_length, chosen_nx, chosen_time, chosen_nt, alpha, function_temperature)

    x, wa = heat_equation_analytical(chosen_length, chosen_nx, chosen_time, chosen_nt, alpha)

    # Plot solutions
    plot_solutions(x, w, wa, chosen_nt, chosen_time)
    plot_surface_solution(x, w, chosen_nt, chosen_time)
