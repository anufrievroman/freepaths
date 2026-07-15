"""Module that controls recording, calculation, and saving thermal and scattering maps"""

from math import cos, sin
from scipy.constants import hbar, pi, k as k_B
import numpy as np
from freepaths.config import cf
from freepaths.options import SimulationMode

class Maps:
    """Parent maps class with functions common to all classes below"""
    def read_data(self, data_dict):
        """Read the data from the finished worker and add new data to already existing"""
        for key, value in data_dict.items():
            setattr(self, key, getattr(self, key) + value)


class ScatteringMap(Maps):
    """Map of scattering in the structure"""

    def __init__(self):
        """Initialize arrays of scattering maps"""
        self.diffuse_scattering_map_x = []
        self.diffuse_scattering_map_y = []
        self.specular_scattering_map_x = []
        self.specular_scattering_map_y = []
        self.internal_scattering_map_x = []
        self.internal_scattering_map_y = []

    def add_scattering_to_map(self, pt, scattering_types):
        """Record the place where a scattering event occurred according to the event type"""

        # Diffuse surface scattering:
        if scattering_types.is_diffuse:
            self.diffuse_scattering_map_x.append(pt.x)
            self.diffuse_scattering_map_y.append(pt.y)

        # Internal scattering:
        elif scattering_types.is_internal:
            self.internal_scattering_map_x.append(pt.x)
            self.internal_scattering_map_y.append(pt.y)

        # Specular surface scattering:
        else:
            self.specular_scattering_map_x.append(pt.x)
            self.specular_scattering_map_y.append(pt.y)

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
        np.savetxt("Data/Scattering map.csv", data, fmt='%1.3e', delimiter=",", header=header, encoding='utf-8')

    def dump_data(self):
        """Return data of a process in the form of a dictionary to be attached to the global data"""
        return {
            'diffuse_scattering_map_x': self.diffuse_scattering_map_x,
            'diffuse_scattering_map_y': self.diffuse_scattering_map_y,
            'specular_scattering_map_x': self.specular_scattering_map_x,
            'specular_scattering_map_y': self.specular_scattering_map_y,
            'internal_scattering_map_x': self.internal_scattering_map_x,
            'internal_scattering_map_y': self.internal_scattering_map_y,
        }


class ThermalMaps(Maps):
    """Maps and profiles of thermal energy in the structure"""

    def __init__(self):
        """Initialize arrays of thermal maps and other parameters"""
        self.thermal_map = np.zeros((cf.number_of_pixels_y, cf.number_of_pixels_x))
        self.effective_heat_flux_profile_y = np.zeros((cf.number_of_pixels_y, cf.number_of_timeframes))
        self.material_heat_flux_profile_y = np.zeros((cf.number_of_pixels_y, cf.number_of_timeframes))
        self.temperature_profile_y = np.zeros((cf.number_of_pixels_y, cf.number_of_timeframes))
        self.heat_flux_map_x = np.zeros((cf.number_of_pixels_y, cf.number_of_pixels_x))
        self.heat_flux_map_y = np.zeros((cf.number_of_pixels_y, cf.number_of_pixels_x))
        self.heat_flux_map_xy = np.zeros((cf.number_of_pixels_y, cf.number_of_pixels_x))
        # Timeframes span the whole virtual time over which phonon emission is spread,
        # so that every phonon contributes to the profiles regardless of its start time:
        self.timesteps_per_timeframe = cf.number_of_virtual_timesteps // cf.number_of_timeframes
        self.effective_thermal_conductivity = np.zeros((cf.number_of_timeframes, 2))
        self.material_thermal_conductivity = np.zeros((cf.number_of_timeframes, 2))

        # Initialize array for thermal conductivity in each time interval
        self.effective_thermal_conductivity[:, 0] = range(cf.number_of_timeframes)
        self.effective_thermal_conductivity[:, 0] += 0.5
        self.effective_thermal_conductivity[:, 0] *= self.timesteps_per_timeframe * cf.timestep * 1e9
        self.material_thermal_conductivity[:, 0] = self.effective_thermal_conductivity[:, 0]

        # Calculate the volumes [m^3] and other parameters (need to be corrected with volume of the holes):
        self.vol_cell_y = cf.length * cf.thickness * cf.width / cf.number_of_pixels_y
        self.vol_pixel =  cf.length * cf.thickness * cf.width / (cf.number_of_pixels_x * cf.number_of_pixels_y)

        # Calculate the pixel volumes with respect to holes:
        self.vol_pixel_ratio = self.calculate_pixel_volumes(cf.number_of_pixels_x, cf.number_of_pixels_y)
        self.vol_pixel_correction_per_row = np.mean(self.vol_pixel_ratio, axis=1)

    def calculate_pixel_volumes(self, number_of_pixels_x, number_of_pixels_y):
        """Calculate a map showing if the pixel contains material (1) or a hole (0)"""
        if not cf.holes:
            return np.ones((number_of_pixels_y, number_of_pixels_x))
        pixel_volume_ratios = np.zeros((number_of_pixels_y, number_of_pixels_x))
        for x_index in range(number_of_pixels_x):
            for y_index in range(number_of_pixels_y):

                # calculate pixel coordinates of pixel center
                y_coord = cf.length / cf.number_of_pixels_y * (y_index + 0.5)
                x_coord = -cf.width/2 + cf.width / cf.number_of_pixels_x * (x_index + 0.5)

                pixel_volume_ratios[y_index, x_index] = not any(hole.is_inside(x_coord, y_coord, None, cf) for hole in cf.holes)

        return pixel_volume_ratios

    def add_energy_to_maps(self, pt, timestep_number, material):
        """
        Register the phonon in the pixel corresponding to its current position at certain timestep
        and adds it to thermal maps and profiles.
        Note that in case of temperature and heat flux profiles, phonon is registered not at
        its current timestep number but at a virtual timestep which is current time + first timestep.
        This first time step is unique for each phonon. This is done to simulate more relistic continious heat flow.
        """

        # Calculate the index of the pixel in which we record this phonon:
        index_x = int(((pt.x + cf.width / 2) * cf.number_of_pixels_x) // cf.width)
        index_y = int(pt.y // (cf.length / cf.number_of_pixels_y))

        # Prevent error if the phonon is outside the structure:
        if (0 <= index_x < cf.number_of_pixels_x) and (0 <= index_y < cf.number_of_pixels_y):

            # Calculate pixel volume correction factors:
            vol_pixel_correction = self.vol_pixel_ratio[index_y, index_x]
            vol_pixel_correction_y = self.vol_pixel_correction_per_row[index_y]

            # Do not record data if the pixel is an empty one:
            if vol_pixel_correction == 0 and cf.ignore_faulty_particles:
                return

            # Record energy [J] and heat flux [W/s/m^2] of this phonon into the pixel of thermal map.
            # With rethermalization, particles are equal-energy bundles (Peraud & Hadjiconstantinou,
            # PRB 84, 205331 (2011)): the deposited weight must not depend on the current frequency,
            # so that the ensemble energy is conserved when frequencies are re-drawn. The absolute
            # value cancels in the thermal conductivity; k_B*T is used as a natural scale.
            if cf.sample_from_dispersion:
                energy = k_B * cf.temp
            else:
                energy = hbar * 2 * pi * pt.f
            self.thermal_map[index_y, index_x] += energy
            self.heat_flux_map_x[index_y, index_x] += energy * sin(pt.theta) * abs(cos(pt.phi)) * pt.speed / self.vol_pixel
            self.heat_flux_map_y[index_y, index_x] += energy * cos(pt.theta) * abs(cos(pt.phi)) * pt.speed / self.vol_pixel

            # Calculate to which timeframe this timestep belongs:
            timeframe_number = (pt.first_timestep + timestep_number) // self.timesteps_per_timeframe

            # Record temperature and heat flux into the corresponding time segment.
            # By default, use the dispersion-only heat capacity, self-consistent with the
            # dispersion-based phonon sampling; the real (experimental) heat capacity,
            # which also counts branches absent from the dispersion, is optional:
            volumetric_heat_capacity = (material.dispersion_heat_capacity if cf.use_dispersion_heat_capacity
                                        else material.heat_capacity * material.density)
            if timeframe_number < cf.number_of_timeframes and vol_pixel_correction_y != 0:
                self.effective_heat_flux_profile_y[index_y, timeframe_number] += energy * cos(pt.theta) * abs(cos(pt.phi)) * pt.speed / self.vol_cell_y
                self.material_heat_flux_profile_y[index_y, timeframe_number] += energy * cos(pt.theta) * abs(cos(pt.phi)) * pt.speed / self.vol_cell_y / vol_pixel_correction_y
                self.temperature_profile_y[index_y, timeframe_number] += energy / volumetric_heat_capacity / self.vol_cell_y / vol_pixel_correction_y


    def calculate_heat_flux_modulus(self):
        """Calculate heat flux modulus as sqrt(q_x^2 + q_y^2)"""
        self.heat_flux_map_xy += np.sqrt(self.heat_flux_map_x**2 + self.heat_flux_map_y**2)


    def calculate_thermal_conductivity(self):
        """Calculate the thermal conductivity for each time interval from heat flux
        and temperature profiles accumulated in that interval"""

        # Restrict the fit to a fraction of the length (GRADIENT_FIT_RANGE) to exclude
        # the quasi-ballistic contact regions near the hot and cold sides, where the
        # temperature profile deviates from linear (temperature jumps at the contacts):
        fit_start = int(cf.gradient_fit_range[0] * cf.number_of_pixels_y)
        fit_end = max(int(cf.gradient_fit_range[1] * cf.number_of_pixels_y), fit_start + 2)

        # For each time interval calculate the thermal conductivity:
        for timeframe_number in range(cf.number_of_timeframes):

            # ATTENTION: This only works when the hot side is at the bottom!
            # We need to improve the d_T calculation so that the gradient is calculated
            # from hot to cold in any direction.

            # Temperature gradient obtained by linear fit of T(y):
            coordinates_y = np.arange(cf.number_of_pixels_y) * cf.length / cf.number_of_pixels_y
            slope, _ = np.polyfit(coordinates_y[fit_start:fit_end],
                                  self.temperature_profile_y[fit_start:fit_end, timeframe_number], 1)
            grad_T = -slope

            # Average heat flux over the same range (always skip first pixel — it has a spurious
            # spike due to rethermalization being applied before map recording in the time-step loop):
            flux_start = max(fit_start, 1)
            J_effective = np.mean(self.effective_heat_flux_profile_y[flux_start:fit_end, timeframe_number])
            J_material = np.mean(self.material_heat_flux_profile_y[flux_start:fit_end, timeframe_number])

            # By definition, J = -K * grad(T), so:
            if grad_T != 0:
                self.effective_thermal_conductivity[timeframe_number, 1] = J_effective / grad_T
                self.material_thermal_conductivity[timeframe_number, 1] = J_material / grad_T

            # Calculate single averaged value in the steady state range:
            self.av_effective_thermal_conductivity = np.mean(self.effective_thermal_conductivity[cf.number_of_stabilization_timeframes:, 1])
            self.av_material_thermal_conductivity = np.mean(self.material_thermal_conductivity[cf.number_of_stabilization_timeframes:, 1])
            self.std_effective_thermal_conductivity = np.std(self.effective_thermal_conductivity[cf.number_of_stabilization_timeframes:, 1])
            self.std_material_thermal_conductivity = np.std(self.material_thermal_conductivity[cf.number_of_stabilization_timeframes:, 1])


    def write_into_files(self, mode=None):
        """Write thermal maps into files"""
        if mode is not SimulationMode.PHONON_TRACING:
            return

        # Create coordinate array [um]
        num_of_points_y = self.temperature_profile_y.shape[0]
        coordinates_y = np.arange(num_of_points_y) * 1e6 * cf.length / num_of_points_y

        # Saving profiles:
        data_temp_y = np.vstack((coordinates_y, self.temperature_profile_y.T)).T
        data_flux_y = np.vstack((coordinates_y, self.effective_heat_flux_profile_y.T, self.material_heat_flux_profile_y.T)).T

        # Saving the thermal conductivity data:
        data_tc = np.vstack((self.effective_thermal_conductivity.T, self.material_thermal_conductivity[:, 1])).T
        av_start = cf.number_of_stabilization_timeframes * self.timesteps_per_timeframe * cf.timestep * 1e9
        av_end = cf.number_of_timeframes * self.timesteps_per_timeframe * cf.timestep * 1e9
        data_av_tc = np.vstack((self.av_effective_thermal_conductivity, self.av_material_thermal_conductivity,
                                self.std_effective_thermal_conductivity, self.std_material_thermal_conductivity, av_start, av_end)).T

        t_headers = ', '.join([f'T (K) [step{i+1}]' for i in range(cf.number_of_timeframes)])
        np.savetxt("Data/Temperature profiles y.csv", data_temp_y, fmt='%1.3e', delimiter=",", header="Y (um), " + t_headers, encoding='utf-8')
        j_headers = ', '.join([f'J_eff (a.u.) [step{i+1}]' for i in range(cf.number_of_timeframes)] + [f'J_mat (a.u.) [step{i+1}]' for i in range(cf.number_of_timeframes)])
        np.savetxt("Data/Heat flux profiles y.csv", data_flux_y, fmt='%1.3e', delimiter=",", header="Y (um), " + j_headers, encoding='utf-8')
        np.savetxt("Data/Thermal conductivity.csv", data_tc, fmt='%1.3e', delimiter=",", header="t(ns), K_eff (W/mK), K_mat (W/mK)", encoding='utf-8')
        np.savetxt("Data/Average thermal conductivity.csv", data_av_tc, fmt='%1.3e', delimiter=",",
                   header="K_eff (W/mK), K_mat (W/mK), error_eff (W/mK), error_mat (W/mK), Av. start (ns), Av. end (ns)", encoding='utf-8')
        np.savetxt("Data/Pixel volumes.csv", self.vol_pixel_ratio, fmt='%1i', delimiter=",", encoding='utf-8')
        np.savetxt("Data/Thermal map.csv", self.thermal_map, fmt='%1.2e', delimiter=",", encoding='utf-8')
        np.savetxt("Data/Heat flux map xy.csv", self.heat_flux_map_xy, fmt='%1.2e', delimiter=",", encoding='utf-8')
        np.savetxt("Data/Heat flux map x.csv", self.heat_flux_map_x, fmt='%1.2e', delimiter=",", encoding='utf-8')
        np.savetxt("Data/Heat flux map y.csv", self.heat_flux_map_y, fmt='%1.2e', delimiter=",", encoding='utf-8')

    def dump_data(self):
        """Return data of a process in the form of a dictionary to be attached to the global data"""
        return {
            'thermal_map': self.thermal_map,
            'heat_flux_map_x': self.heat_flux_map_x,
            'heat_flux_map_y': self.heat_flux_map_y,
            'heat_flux_map_xy': self.heat_flux_map_xy,
            'effective_heat_flux_profile_y': self.effective_heat_flux_profile_y,
            'material_heat_flux_profile_y': self.material_heat_flux_profile_y,
            'temperature_profile_y': self.temperature_profile_y,
        }