import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Provided data
temperature = np.array([4.0, 21.5, 47.7, 69.2, 100.0, 123.3, 150.0, 175.1, 200.0, 224.5, 273.0, 293.7, 300.0, 315.2, 335.1, 395.6, 413.9, 423.5, 455.3, 475.9, 494.2, 515.7, 534.8, 574.5])
heat_capacity = np.array([0.65, 2.5, 10.1, 25.4, 71.3, 109.6, 228.3, 349.6, 415.0, 464.4, 566.5, 618.8, 627.0, 645.6, 700.4, 808.9, 830.5, 845.8, 885.4, 900.6, 924.9, 941.4, 950.3, 969.4])

# Splitting data into three parts
below_90K_indices = temperature < 90
between_90_and_200K_indices = (temperature >= 90) & (temperature <= 200)
above_200K_indices = temperature > 200

below_90K_temperature = temperature[below_90K_indices]
between_90_and_200K_temperature = temperature[between_90_and_200K_indices]
above_200K_temperature = temperature[above_200K_indices]

below_90K_heat_capacity = heat_capacity[below_90K_indices]
between_90_and_200K_heat_capacity = heat_capacity[between_90_and_200K_indices]
above_200K_heat_capacity = heat_capacity[above_200K_indices]

# Define the functions to fit (3rd-degree polynomial for each part)
def third_degree_func(x, a, b, c, d):
    return a * x**3 + b * x**2 + c * x + d

# Perform curve fitting for each part
below_90K_popt, _ = curve_fit(third_degree_func, below_90K_temperature, below_90K_heat_capacity)
between_90_and_200K_popt, _ = curve_fit(third_degree_func, between_90_and_200K_temperature, between_90_and_200K_heat_capacity)
above_200K_popt, _ = curve_fit(third_degree_func, above_200K_temperature, above_200K_heat_capacity)

# Generate y-values for the fitted curves
below_90K_y_fit = third_degree_func(below_90K_temperature, *below_90K_popt)
between_90_and_200K_y_fit = third_degree_func(between_90_and_200K_temperature, *between_90_and_200K_popt)
above_200K_y_fit = third_degree_func(above_200K_temperature, *above_200K_popt)

# Plotting
plt.figure(figsize=(10, 6))

# Plot below 90K fit
plt.scatter(below_90K_temperature, below_90K_heat_capacity, label='Below 90 K', color='blue')
plt.plot(below_90K_temperature, below_90K_y_fit, color='blue', linestyle='--')

# Plot between 90 and 200K fit
plt.scatter(between_90_and_200K_temperature, between_90_and_200K_heat_capacity, label='Between 90 K and 200 K', color='orange')
plt.plot(between_90_and_200K_temperature, between_90_and_200K_y_fit, color='orange', linestyle='--')

# Plot above 200K fit
plt.scatter(above_200K_temperature, above_200K_heat_capacity, label='Above 200 K', color='green')
plt.plot(above_200K_temperature, above_200K_y_fit, color='green', linestyle='--')

plt.xlabel('Temperature (K)')
plt.ylabel('Heat Capacity')
plt.title('Fitted Curve of Heat Capacity vs Temperature (Log-Log Scale)')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.grid(True)

plt.show()

# Show the fitted parameters
print("Fitted parameters for below 90 K (a, b, c, d):", below_90K_popt)
print("Fitted parameters for between 90 K and 200 K (a, b, c, d):", between_90_and_200K_popt)
print("Fitted parameters for above 200 K (a, b, c, d):", above_200K_popt)

# Function to calculate heat capacity using fitting coefficients
def calculate_heat_capacity(temperature):
    if temperature < 90:
        coeffs = below_90K_popt
    elif 90 <= temperature <= 200:
        coeffs = between_90_and_200K_popt
    else:
        coeffs = above_200K_popt

    # Third-degree polynomial calculation
    return np.polyval(coeffs, temperature)

# Example usage
temperature = 150  # Example temperature
heat_capacity = calculate_heat_capacity(temperature)
print(f"Heat capacity at {temperature} K:", heat_capacity)
