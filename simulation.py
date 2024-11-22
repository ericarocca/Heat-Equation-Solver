import configparser
import numpy as np
from function import heat_equation_CN, heat_equation_analytical, function_temperature, plot_solutions, plot_surface_solution, check_stability

def process_configuration(config_file):
    """Function that process a single configuration file."""
    config = configparser.ConfigParser()
    config.read(config_file)

    length = float(config.get('settings', 'length'))
    nx_values = list(map(int, config.get('settings', 'nx_values').split(',')))
    time = float(config.get('settings', 'time'))
    nt_values = list(map(int, config.get('settings', 'nt_values').split(',')))
    alpha = float(config.get('settings', 'alpha'))

    numerical_solution = config.get('paths', 'numerical_solution')
    analytical_solution = config.get('paths', 'analytical_solution')

    #verify the presence of stable combinations, then solve and plot for those
    stable_combinations = check_stability(length, time, alpha, nx_values, nt_values)
    if not stable_combinations:
        raise ValueError(f"No stable combinations found for parameters in {config_file}.")

    for combination in stable_combinations:
        chosen_length, chosen_time, chosen_nx, chosen_nt, chosen_r = combination

        print(f"Running simulation with nx={chosen_nx}, nt={chosen_nt}, r={chosen_r}")

        x, w = heat_equation_CN(chosen_length, chosen_nx, chosen_time, chosen_nt, alpha, function_temperature)
        x, wa = heat_equation_analytical(chosen_length, chosen_nx, chosen_time, chosen_nt, alpha)

        np.save(numerical_solution, w)
        np.save(analytical_solution, wa)

        plot_solutions(x, w, wa, chosen_nt, chosen_time, chosen_length, chosen_nx, alpha)
        plot_surface_solution(x, w, chosen_nt, chosen_time, chosen_length, chosen_nx, alpha)

#the user can choose a specific configuration
print("Choose configuration A or B:")
choice = input("Enter A or B: ").strip().upper()

#map choices to configuration files
configuration_map = {
    "A": "configurationA.txt",
    "B": "configurationB.txt"
}

#check if the choice is valid
if choice not in configuration_map:
    raise ValueError("Invalid choice!!! Please enter 'A' or 'B'.")

#process the chosen configuration
process_configuration(configuration_map[choice])
