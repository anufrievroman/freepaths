"""Module that controles recording and calculation of maps"""

import random
from math import cos
from scipy.constants import hbar, pi
import numpy as np

from freepaths.config import cf


class ScatteringMap:
    """Map of scattering in the structure"""

    def __init__(self):
        """Initialize arrays of scattering maps"""
        self.diffuse_scattering_map_x = []
        self.diffuse_scattering_map_y = []
        self.specular_scattering_map_x = []
        self.specular_scattering_map_y = []
        self.internal_scattering_map_x = []
        self.internal_scattering_map_y = []

    def add_scattering_to_map(self, ph, scattering_types):
        """Record the place where a scattering event occurred according to the event type"""

        # Diffuse surface scattering:
        if scattering_types.is_diffuse:
            self.diffuse_scattering_map_x.append(ph.x)
            self.diffuse_scattering_map_y.append(ph.y)

        # Internal scattering:
        elif scattering_types.is_internal:
            self.internal_scattering_map_x.append(ph.x)
            self.internal_scattering_map_y.append(ph.y)

        # Specular surface scattering:
        else:
            self.specular_scattering_map_x.append(ph.x)
            self.specular_scattering_map_y.append(ph.y)

    def write_into_files(self):
        """Write scattering map into file"""

        # Maximal number of scattering event of any type:
        n_max = max([len(self.specular_scattering_map_x),
                     len(self.internal_scattering_map_x),
                     len(self.diffuse_scattering_map_x)])

        # Create an array and fill it with the coordinates:
        data = np.zeros((n_max, 6))
        for index, (x, y) in enumerate(zip(self.specular_scattering_map_x, self.specular_scattering_map_y)):
            data[index, 0] = x
            data[index, 1] = y
        for index, (x, y) in enumerate(zip(self.diffuse_scattering_map_x, self.diffuse_scattering_map_y)):
            data[index, 2] = x
            data[index, 3] = y
        for index, (x, y) in enumerate(zip(self.internal_scattering_map_x, self.internal_scattering_map_y)):
            data[index, 4] = x
            data[index, 5] = y

        # Save into file:
        header = "Specular X, Specular Y, Diffuse X, Diffuse Y, Internal X, Internal Y"
        np.savetxt("Data/Scattering map.csv", data, fmt='%1.2e', delimiter=",", header=header, encoding='utf-8')


class ThermalMaps:
    """Maps and profiles of thermal energy in the structure"""

    def __init__(self):
        """Initialize arrays of thermal maps"""
        self.thermal_map = np.zeros((cf.number_of_pixels_y, cf.number_of_pixels_x))
        self.heat_flux_profile_x = np.zeros((cf.number_of_pixels_x, cf.number_of_timeframes))
        self.heat_flux_profile_y = np.zeros((cf.number_of_pixels_y, cf.number_of_timeframes))
        self.temperature_profile_x = np.zeros((cf.number_of_pixels_x, cf.number_of_timeframes))
        self.temperature_profile_y = np.zeros((cf.number_of_pixels_y, cf.number_of_timeframes))
        self.thermal_conductivity = np.zeros((cf.number_of_timeframes, 2))

    def add_energy_to_maps(self, ph, timestep_number, material):
        """This function registers the phonon in the pixel corresponding to its current position
        and at certain timesteps and adds it to thermal maps and thermal profiles"""

        # Calculate the index of the pixel in which this phonon is now:
        index_x = int(((ph.x + cf.width / 2) * cf.number_of_pixels_x) // cf.width)
        # index_y = int((ph.y*number_of_pixels_y) // length)
        index_y = int(ph.y // (cf.length / cf.number_of_pixels_y))

        # Calculate the volume of this pixel:
        vol_cell = cf.length * cf.thickness * cf.width
        vol_cell_x = vol_cell / cf.number_of_pixels_x
        vol_cell_y = vol_cell / cf.number_of_pixels_y

        # Here we arbitrarily correct the volume of the unit cells in pillars:
        if cf.include_pillars == 'yes':
            vol_cell_x += 2.5 * 0.3333 * cf.pillar_height * (cf.circular_hole_diameter / 2) ** 2
            vol_cell_y += 2.5 * 0.3333 * cf.pillar_height * (cf.circular_hole_diameter / 2) ** 2

        # Prevent error if the phonon is outside the structure:
        if (0 <= index_x < cf.number_of_pixels_x) and (0 <= index_y < cf.number_of_pixels_y):

            # Record energy h*w of this phonon into the pixel of thermal map:
            energy = hbar * 2 * pi * ph.f
            self.thermal_map[index_y, index_x] += energy

            # Record energy of this phonon into flux and temperature profiles: (DOUBLE-CHECK THIS)
            random_timeframe = random.randint(0, cf.number_of_timesteps)
            assigned_time = (timestep_number + random_timeframe) * cf.timestep * cf.number_of_timeframes
            total_time = cf.number_of_timesteps * cf.timestep
            timeframe_number = int(assigned_time // total_time)

            if timeframe_number < cf.number_of_timeframes:
                self.heat_flux_profile_x[index_x, timeframe_number] += energy * cos(ph.theta) * abs(cos(ph.phi)) * ph.speed / vol_cell_x
                self.heat_flux_profile_y[index_y, timeframe_number] += energy * cos(ph.theta) * abs(cos(ph.phi)) * ph.speed / vol_cell_y
                self.temperature_profile_x[index_x, timeframe_number] += energy / (cf.specific_heat_capacity * material.density) / vol_cell_x
                self.temperature_profile_y[index_y, timeframe_number] += energy / (cf.specific_heat_capacity * material.density) / vol_cell_y

    def calculate_thermal_conductivity(self):
        """Calculate the thermal conductivity for each time interval from heat flux
        and temperature profiles accumulated in that interval"""

        # Initialize array for thermal conductivity in each time interval
        self.thermal_conductivity[:, 0] = range(cf.number_of_timeframes)
        total_time = cf.number_of_timesteps * cf.timestep * 1e9
        self.thermal_conductivity[:, 0] *= total_time / cf.number_of_timeframes

        # For each time interval calculate the thermal conductivity:
        for timeframe_number in range(cf.number_of_timeframes):
            # Here we ignore the first pixel, because there are anomalies usually...
            T_high = self.temperature_profile_y[1, timeframe_number]
            T_low = self.temperature_profile_y[(cf.number_of_pixels_y - 1), timeframe_number]

            # Temperature gradient:
            d_T = T_high - T_low

            # Average heat flux:
            J = sum(self.heat_flux_profile_y[1:cf.number_of_pixels_y, timeframe_number]) / (cf.number_of_pixels_y - 1)

            # Here dL is shorter than actual length because we ignore 1st pixel and lose one more due to averaging:
            d_L = (cf.number_of_pixels_y - 2) * cf.length / cf.number_of_pixels_y

            # By definition, J = -K*grad(T), so the thermal conductivity:
            self.thermal_conductivity[timeframe_number, 1] = J * d_L / d_T

    def write_into_files(self):
        """Write thermal map into file"""

        if cf.output_raw_thermal_map:
            np.savetxt("Data/Thermal map.csv", self.thermal_map, fmt='%1.2e', delimiter=",", encoding='utf-8')

        # Create coordinate arrays [um]
        num_of_points_x = self.temperature_profile_x.shape[0]
        num_of_points_y = self.temperature_profile_y.shape[0]
        coordinates_x = np.arange(num_of_points_x) * 1e6 * cf.width / num_of_points_x
        coordinates_y = np.arange(num_of_points_y) * 1e6 * cf.length / num_of_points_y

        # Saving all the profiles in the files:
        data_temp_x = np.vstack((coordinates_x, self.temperature_profile_x.T)).T
        data_temp_y = np.vstack((coordinates_y, self.temperature_profile_y.T)).T
        data_flux_x = np.vstack((coordinates_x, self.heat_flux_profile_x.T)).T
        data_flux_y = np.vstack((coordinates_y, self.heat_flux_profile_y.T)).T
        data_tc = self.thermal_conductivity
        np.savetxt("Data/Temperature profiles x.csv", data_temp_x, fmt='%1.3e', delimiter=",", header="X (um), T (K)", encoding='utf-8')
        np.savetxt("Data/Temperature profiles y.csv", data_temp_y, fmt='%1.3e', delimiter=",", header="Y (um), T (K)", encoding='utf-8')
        np.savetxt("Data/Heat flux profiles x.csv", data_flux_x, fmt='%1.3e', delimiter=",", header="Y (um), J (a.u.)", encoding='utf-8')
        np.savetxt("Data/Heat flux profiles y.csv", data_flux_y, fmt='%1.3e', delimiter=",", header="Y (um), J (a.u.)", encoding='utf-8')
        np.savetxt("Data/Thermal conductivity.csv", data_tc, fmt='%1.3e', delimiter=",", header="t(ns), K (W/mK)", encoding='utf-8')
