import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Provided data
temperature = np.array([3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 75, 80, 90, 100, 110, 120, 125, 130, 140, 150, 160, 170, 175, 180, 190, 200, 210, 220, 225, 230, 240, 250, 260, 270, 280, 290, 300])
heat_capacity = np.array([0.0074, 0.0176, 0.0344, 0.0594, 0.0944, 0.1410, 0.2007, 0.2756, 1.0874, 3.3669, 8.4921, 17.1266, 28.8410, 44.0449, 60.8866, 78.5117, 115.0080, 151.3619, 169.6635, 188.1075, 224.0698, 259.4267, 294.9973, 329.3573, 346.1278, 362.684, 395.478, 425.672, 454.904, 482.677, 495.673, 508.670, 533.736, 556.952, 577.639, 596.795, 606.622, 615.666, 633.006, 648.709, 663.165, 677.051, 690.083, 702.795, 714.473])

# Splitting data into three parts
below_20K_indices = temperature < 20
between_20_and_50K_indices = (temperature >= 20) & (temperature <= 50)
above_50K_indices = temperature > 50

below_20K_temperature = temperature[below_20K_indices]
between_20_and_50K_temperature = temperature[between_20_and_50K_indices]
above_50K_temperature = temperature[above_50K_indices]

below_20K_heat_capacity = heat_capacity[below_20K_indices]
between_20_and_50K_heat_capacity = heat_capacity[between_20_and_50K_indices]
above_50K_heat_capacity = heat_capacity[above_50K_indices]

# Define the functions to fit (3rd-degree polynomial for each part)
def third_degree_func(x, a, b, c, d):
    return a * x**3 + b * x**2 + c * x + d

# Perform curve fitting for each part
below_20K_popt, _ = curve_fit(third_degree_func, below_20K_temperature, below_20K_heat_capacity)
between_20_and_50K_popt, _ = curve_fit(third_degree_func, between_20_and_50K_temperature, between_20_and_50K_heat_capacity)
above_50K_popt, _ = curve_fit(third_degree_func, above_50K_temperature, above_50K_heat_capacity)

# Generate y-values for the fitted curves
below_20K_y_fit = third_degree_func(below_20K_temperature, *below_20K_popt)
between_20_and_50K_y_fit = third_degree_func(between_20_and_50K_temperature, *between_20_and_50K_popt)
above_50K_y_fit = third_degree_func(above_50K_temperature, *above_50K_popt)

# Plotting
plt.figure(figsize=(10, 6))

# Plot below 20K fit
plt.scatter(below_20K_temperature, below_20K_heat_capacity, label='Below 20 K', color='blue')
plt.plot(below_20K_temperature, below_20K_y_fit, color='blue', linestyle='--')

# Plot between 20 and 50K fit
plt.scatter(between_20_and_50K_temperature, between_20_and_50K_heat_capacity, label='Between 20 K and 50 K', color='orange')
plt.plot(between_20_and_50K_temperature, between_20_and_50K_y_fit, color='orange', linestyle='--')

# Plot above 50K fit
plt.scatter(above_50K_temperature, above_50K_heat_capacity, label='Above 50 K', color='green')
plt.plot(above_50K_temperature, above_50K_y_fit, color='green', linestyle='--')

plt.xlabel('Temperature (K)')
plt.ylabel('Heat Capacity (J/kg/K)')
plt.title('Fitted Curve of Heat Capacity vs Temperature')
plt.xscale('log')
plt.yscale('log')
plt.legend()

plt.show()

# Show the fitted parameters
print("Fitted parameters for below 20 K (a, b, c, d):", below_20K_popt)
print("Fitted parameters for between 20 K and 50 K (a, b, c, d):", between_20_and_50K_popt)
print("Fitted parameters for above 50 K (a, b, c, d):", above_50K_popt)
