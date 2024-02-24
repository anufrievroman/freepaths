import numpy as np
import matplotlib.pyplot as plt

# Define the equation
def heat_capacity(temperature):
    coeffs = np.array([6.309e-9, 6.27e-6, 8.729e-4, 0])
    return 1000 * np.polyval(coeffs, temperature)

# Generate temperature values from 4 K to 300 K
temperature_range = np.linspace(80, 300, 1000)

# Calculate heat capacity for the temperature range
heat_capacity_values = heat_capacity(temperature_range)

# Plot the equation
plt.plot(temperature_range, heat_capacity_values, label='Heat Capacity')
plt.xlabel('Temperature (K)')
plt.ylabel('Heat Capacity (J/kg/K)')
plt.title('Heat Capacity vs Temperature')
plt.grid(True)
plt.legend()
plt.show()
