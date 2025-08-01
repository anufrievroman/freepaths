"""Module that controls recording, calculation, and saving thermal and scattering maps"""

from math import cos
from scipy.constants import hbar, pi, elementary_charge
import numpy as np
from math import cos , sin
from freepaths.config import cf
from freepaths.move import move

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
        self.effective_heat_flux_profile_x = np.zeros((cf.number_of_pixels_x, cf.number_of_timeframes))
        self.effective_heat_flux_profile_y = np.zeros((cf.number_of_pixels_y, cf.number_of_timeframes))
        self.material_heat_flux_profile_x = np.zeros((cf.number_of_pixels_x, cf.number_of_timeframes))
        self.material_heat_flux_profile_y = np.zeros((cf.number_of_pixels_y, cf.number_of_timeframes))
        self.temperature_profile_x = np.zeros((cf.number_of_pixels_x, cf.number_of_timeframes))
        self.temperature_profile_y = np.zeros((cf.number_of_pixels_y, cf.number_of_timeframes))
        self.heat_flux_map_x = np.zeros((cf.number_of_pixels_y, cf.number_of_pixels_x))
        self.heat_flux_map_y = np.zeros((cf.number_of_pixels_y, cf.number_of_pixels_x))
        self.heat_flux_map_xy = np.zeros((cf.number_of_pixels_y, cf.number_of_pixels_x))
        self.heat_flux_map_x_weighted = np.zeros((cf.number_of_pixels_y, cf.number_of_pixels_x))
        self.heat_flux_map_y_weighted = np.zeros((cf.number_of_pixels_y, cf.number_of_pixels_x))
        self.number_phonons_in_pixel = np.zeros((cf.number_of_pixels_y, cf.number_of_pixels_x))
        self.timesteps_per_timeframe = cf.number_of_timesteps // cf.number_of_timeframes
        self.effective_thermal_conductivity = np.zeros((cf.number_of_timeframes, 2))
        self.material_thermal_conductivity = np.zeros((cf.number_of_timeframes, 2))

        # Initialize array for thermal conductivity in each time interval
        self.effective_thermal_conductivity[:, 0] = range(cf.number_of_timeframes)
        self.effective_thermal_conductivity[:, 0] += 0.5
        self.effective_thermal_conductivity[:, 0] *= self.timesteps_per_timeframe * cf.timestep * 1e9
        self.material_thermal_conductivity[:, 0] = self.effective_thermal_conductivity[:, 0]

        # Calculate the volumes [m^3] and other parameters (need to be corrected with volume of the holes):
        self.vol_cell_x = cf.length * cf.thickness * cf.width / cf.number_of_pixels_x
        self.vol_cell_y = cf.length * cf.thickness * cf.width / cf.number_of_pixels_y
        self.vol_pixel =  cf.length * cf.thickness * cf.width / (cf.number_of_pixels_x * cf.number_of_pixels_y)

        # Calculate the pixel volumes with respect to holes:
        self.vol_pixel_ratio = self.calculate_pixel_volumes(cf.number_of_pixels_x, cf.number_of_pixels_y)

    def calculate_pixel_volumes(self, number_of_pixels_x, number_of_pixels_y):
        """Calculate a map showing if the pixel contains material (1) or a hole (0)"""
        pixel_volume_ratios = np.zeros((number_of_pixels_y, number_of_pixels_x))
        for x_index in range(number_of_pixels_x):
            for y_index in range(number_of_pixels_y):

                # calculate pixel coordinates of pixel center
                y_coord = cf.length / cf.number_of_pixels_y * (y_index + 0.5)
                x_coord = -cf.width/2 + cf.width / cf.number_of_pixels_x * (x_index + 0.5)

                pixel_volume_ratios[y_index, x_index] = not any(hole.is_inside(x_coord, y_coord, None, cf) for hole in cf.holes)

        return pixel_volume_ratios

    def add_energy_to_maps(self, ph, timestep_number, material):
        """
        Register the phonon in the pixel corresponding to its current position at certain timestep
        and adds it to thermal maps and profiles.
        Note that in case of temperature and heat flux profiles, phonon is registered not at
        its current timestep number but at a virtual timestep which is current time + first timestep.
        This first time step is unique for each phonon. This is done to simulate more relistic continious heat flow.
        """

        # Move the recording point forwards half a timestep as to remove asymetric boundary fluxes (doesn't work)
        ph_x, ph_y, _ = move(ph, cf.timestep/2)

        # Calculate the index of the pixel in which we record this phonon:
        index_x = int(((ph_x + cf.width / 2) * cf.number_of_pixels_x) // cf.width)
        index_y = int(ph_y // (cf.length / cf.number_of_pixels_y))

        # Prevent error if the phonon is outside the structure:
        if (0 <= index_x < cf.number_of_pixels_x) and (0 <= index_y < cf.number_of_pixels_y):

            # Calculate pixel volume correction factors:
            vol_pixel_correction = self.vol_pixel_ratio[index_y, index_x]
            vol_pixel_correction_x = np.mean(self.vol_pixel_ratio[:, index_x])
            vol_pixel_correction_y = np.mean(self.vol_pixel_ratio[index_y, :])

            # Do not record data if the pixel is an empty one:
            if vol_pixel_correction == 0 and cf.ignore_faulty_particles:
                return

            # Record energy h*w [J] and heat flux [W/s/m^2] of this phonon into the pixel of thermal map:
            energy = hbar * 2 * pi * ph.f
            self.number_phonons_in_pixel[index_y, index_x] += 1
            self.thermal_map[index_y, index_x] += energy
            self.heat_flux_map_x[index_y, index_x] += energy * sin(ph.theta) * abs(cos(ph.phi)) * ph.speed / self.vol_pixel
            self.heat_flux_map_y[index_y, index_x] += energy * cos(ph.theta) * abs(cos(ph.phi)) * ph.speed / self.vol_pixel

            # Calculate to which timeframe this timestep belongs:
            timeframe_number = (ph.first_timestep + timestep_number) // self.timesteps_per_timeframe

            # Record temperature and energy into the corresponding time segment:
            if timeframe_number < cf.number_of_timeframes and vol_pixel_correction_x != 0 and vol_pixel_correction_y != 0:
                self.effective_heat_flux_profile_x[index_x, timeframe_number] += energy * sin(ph.theta) * abs(cos(ph.phi)) * ph.speed / self.vol_cell_x
                self.effective_heat_flux_profile_y[index_y, timeframe_number] += energy * cos(ph.theta) * abs(cos(ph.phi)) * ph.speed / self.vol_cell_y
                # the material heat_flux_profile could also be calculated afterwards by dividing the effective profile with the volume ratio
                self.material_heat_flux_profile_x[index_x, timeframe_number] += energy * sin(ph.theta) * abs(cos(ph.phi)) * ph.speed / self.vol_cell_x / vol_pixel_correction_x
                self.material_heat_flux_profile_y[index_y, timeframe_number] += energy * cos(ph.theta) * abs(cos(ph.phi)) * ph.speed / self.vol_cell_y / vol_pixel_correction_y

                self.temperature_profile_x[index_x, timeframe_number] += energy / (material.heat_capacity * material.density) / self.vol_cell_x / vol_pixel_correction_x
                self.temperature_profile_y[index_y, timeframe_number] += energy / (material.heat_capacity * material.density) / self.vol_cell_y / vol_pixel_correction_y


    def calculate_weighted_flux(self):
        """Calculate heat flux normalized by the number of phonons in each pixel, except where the number is zero"""
        zero_mask = self.number_phonons_in_pixel == 0
        self.heat_flux_map_x_weighted = np.divide(self.heat_flux_map_x, self.number_phonons_in_pixel,
                                              out=np.zeros_like(self.heat_flux_map_x), where=~zero_mask)
        self.heat_flux_map_y_weighted = np.divide(self.heat_flux_map_y, self.number_phonons_in_pixel,
                                              out=np.zeros_like(self.heat_flux_map_y), where=~zero_mask)


    def calculate_heat_flux_modulus(self):
        """Calculate heat flux modulus as sqrt(q_x^2 + q_y^2)"""
        self.heat_flux_map_xy += np.sqrt(self.heat_flux_map_x**2 + self.heat_flux_map_y**2)


    def calculate_thermal_conductivity(self):
        """Calculate the thermal conductivity for each time interval from heat flux
        and temperature profiles accumulated in that interval"""

        # For each time interval calculate the thermal conductivity:
        for timeframe_number in range(cf.number_of_timeframes):

            # ATTENTION: This only works when the hot side is at the bottom!
            # We need to improve the d_T calculation so that the gradient is calculated
            # from hot to cold in any direction.

            # Temperature gradient obtained by linear fit of T(y):
            coordinates_y = np.arange(cf.number_of_pixels_y) * cf.length / cf.number_of_pixels_y
            slope, _ = np.polyfit(coordinates_y, self.temperature_profile_y[:, timeframe_number], 1)
            grad_T = -slope

            # Average heat flux:
            J_effective = np.mean(self.effective_heat_flux_profile_y[:, timeframe_number])
            J_material = np.mean(self.material_heat_flux_profile_y[:, timeframe_number])

            # By definition, J = -K * grad(T), so:
            if grad_T != 0:
                self.effective_thermal_conductivity[timeframe_number, 1] = J_effective / grad_T
                self.material_thermal_conductivity[timeframe_number, 1] = J_material / grad_T

            # Calculate single averaged value in the steady state range:
            self.av_effective_thermal_conductivity = np.mean(self.effective_thermal_conductivity[cf.number_of_stabilization_timeframes:, 1])
            self.av_material_thermal_conductivity = np.mean(self.material_thermal_conductivity[cf.number_of_stabilization_timeframes:, 1])
            self.std_effective_thermal_conductivity = np.std(self.effective_thermal_conductivity[cf.number_of_stabilization_timeframes:, 1])
            self.std_material_thermal_conductivity = np.std(self.material_thermal_conductivity[cf.number_of_stabilization_timeframes:, 1])



    def write_into_files(self):
        """Write thermal maps into files"""

        # Create coordinate arrays [um]
        num_of_points_x = self.temperature_profile_x.shape[0]
        num_of_points_y = self.temperature_profile_y.shape[0]
        coordinates_x = np.arange(num_of_points_x) * 1e6 * cf.width / num_of_points_x
        coordinates_y = np.arange(num_of_points_y) * 1e6 * cf.length / num_of_points_y

        # Saving all the profiles:
        data_temp_x = np.vstack((coordinates_x, self.temperature_profile_x.T)).T
        data_temp_y = np.vstack((coordinates_y, self.temperature_profile_y.T)).T
        data_flux_x = np.vstack((coordinates_x, self.effective_heat_flux_profile_x.T, self.material_heat_flux_profile_x.T)).T
        data_flux_y = np.vstack((coordinates_y, self.effective_heat_flux_profile_y.T, self.material_heat_flux_profile_y.T)).T

        # Saving the thermal conductivity data:
        data_tc = np.vstack((self.effective_thermal_conductivity.T, self.material_thermal_conductivity[:, 1])).T
        av_start = cf.number_of_stabilization_timeframes * self.timesteps_per_timeframe * cf.timestep * 1e9
        av_end = cf.number_of_timeframes * self.timesteps_per_timeframe * cf.timestep * 1e9
        data_av_tc = np.vstack((self.av_effective_thermal_conductivity, self.av_material_thermal_conductivity,
                                self.std_effective_thermal_conductivity, self.std_material_thermal_conductivity, av_start, av_end)).T

        t_headers = ', '.join([f'T (K) [step{i+1}]' for i in range(cf.number_of_timeframes)])
        np.savetxt("Data/Temperature profiles x.csv", data_temp_x, fmt='%1.3e', delimiter=",", header="X (um), " + t_headers, encoding='utf-8')
        np.savetxt("Data/Temperature profiles y.csv", data_temp_y, fmt='%1.3e', delimiter=",", header="Y (um), " + t_headers, encoding='utf-8')
        j_headers = ', '.join([f'J_eff (a.u.) [step{i+1}]' for i in range(cf.number_of_timeframes)] + [f'J_mat (a.u.) [step{i+1}]' for i in range(cf.number_of_timeframes)])
        np.savetxt("Data/Heat flux profiles x.csv", data_flux_x, fmt='%1.3e', delimiter=",", header="Y (um), " + j_headers, encoding='utf-8')
        np.savetxt("Data/Heat flux profiles y.csv", data_flux_y, fmt='%1.3e', delimiter=",", header="Y (um), " + j_headers, encoding='utf-8')
        np.savetxt("Data/Thermal conductivity.csv", data_tc, fmt='%1.3e', delimiter=",", header="t(ns), K_eff (W/mK), K_mat (W/mK)", encoding='utf-8')
        np.savetxt("Data/Average thermal conductivity.csv", data_av_tc, fmt='%1.3e', delimiter=",",
                   header="K_eff (W/mK), K_mat (W/mK), error_eff (W/mK), error_mat (W/mK), Av. start (ns), Av. end (ns)", encoding='utf-8')

        # Saving thermal maps:
        np.savetxt("Data/Pixel volumes.csv", self.vol_pixel_ratio, fmt='%1i', delimiter=",", encoding='utf-8')
        np.savetxt("Data/Thermal map.csv", self.thermal_map, fmt='%1.2e', delimiter=",", encoding='utf-8')
        np.savetxt("Data/Heat flux map xy.csv", self.heat_flux_map_xy, fmt='%1.2e', delimiter=",", encoding='utf-8')
        np.savetxt("Data/Heat flux map x.csv", self.heat_flux_map_x, fmt='%1.2e', delimiter=",", encoding='utf-8')
        np.savetxt("Data/Heat flux map y.csv", self.heat_flux_map_y, fmt='%1.2e', delimiter=",", encoding='utf-8')
        np.savetxt("Data/Heat flux map x weighted.csv", self.heat_flux_map_x_weighted, fmt='%1.2e', delimiter=",", encoding='utf-8')
        np.savetxt("Data/Heat flux map y weighted.csv", self.heat_flux_map_y_weighted, fmt='%1.2e', delimiter=",", encoding='utf-8')

    def dump_data(self):
        """Return data of a process in the form of a dictionary to be attached to the global data"""
        return {
            'thermal_map': self.thermal_map,
            'heat_flux_map_x': self.heat_flux_map_x,
            'heat_flux_map_y': self.heat_flux_map_y,
            'heat_flux_map_xy': self.heat_flux_map_xy,
            'heat_flux_map_y_weighted': self.heat_flux_map_y_weighted,
            'heat_flux_map_x_weighted': self.heat_flux_map_x_weighted,
            'number_phonons_in_pixel': self.number_phonons_in_pixel,
            'effective_heat_flux_profile_x': self.effective_heat_flux_profile_x,
            'effective_heat_flux_profile_y': self.effective_heat_flux_profile_y,
            'material_heat_flux_profile_x': self.material_heat_flux_profile_x,
            'material_heat_flux_profile_y': self.material_heat_flux_profile_y,
            'temperature_profile_x': self.temperature_profile_x,
            'temperature_profile_y': self.temperature_profile_y,
        }

class ElectricalMaps(Maps):
    """Maps and profiles of current density and conductivity via carrier injection"""

    def __init__(self):
        self.current_density_map_x = np.zeros((cf.number_of_pixels_y, cf.number_of_pixels_x))
        self.current_density_map_y = np.zeros((cf.number_of_pixels_y, cf.number_of_pixels_x))
        self.current_density_profile_x = np.zeros((cf.number_of_pixels_x, cf.number_of_timeframes))
        self.current_density_profile_y = np.zeros((cf.number_of_pixels_y, cf.number_of_timeframes))
        self.density_profile_y = np.zeros((cf.number_of_pixels_y, cf.number_of_timeframes))  # [1/m³]
        self.timesteps_per_timeframe = cf.number_of_timesteps // cf.number_of_timeframes
        self.vol_pixel_ratio = self.calculate_pixel_volumes(cf.number_of_pixels_x, cf.number_of_pixels_y)
        self.electron_count_map = np.zeros((cf.number_of_pixels_y, cf.number_of_pixels_x))
        self.current_density_magnitude = np.zeros((cf.number_of_pixels_y, cf.number_of_pixels_x))
        self.conductivity_from_gradient = np.zeros((cf.number_of_timeframes, 2))
        self.avg_sigma_from_gradient = 0
        self.std_sigma_from_gradient = 0

    def calculate_pixel_volumes(self, number_of_pixels_x, number_of_pixels_y):
        pixel_volume_ratios = np.zeros((number_of_pixels_y, number_of_pixels_x))
        for x in range(number_of_pixels_x):
            for y in range(number_of_pixels_y):
                y_coord = cf.length / number_of_pixels_y * (y + 0.5)
                x_coord = -cf.width / 2 + cf.width / number_of_pixels_x * (x + 0.5)
                pixel_volume_ratios[y, x] = not any(hole.is_inside(x_coord, y_coord, None, cf) for hole in cf.holes)
        return pixel_volume_ratios

    def add_electron_to_maps(self, electron, timestep_number, material):
        if electron.vx is None or electron.vy is None:
            return
        
        # Move midpoint to reduce border bias
        x_mid, y_mid, _ = move(electron, cf.timestep/2)

        Nx, Ny = cf.number_of_pixels_x, cf.number_of_pixels_y
        ix = int(((x_mid + cf.width/2) * Nx) // cf.width)
        iy = int(y_mid // (cf.length / Ny))
        if not (0 <= ix < Nx and 0 <= iy < Ny):
            return
        if self.vol_pixel_ratio[iy, ix] == 0 and cf.ignore_faulty_particles:
            return

        # Charge and pixel volume
        q = -elementary_charge
        Vpix = (cf.length * cf.width * cf.thickness) / (Nx * Ny)

        # Velocity components (assumes vx/vy exist)
        vx, vy = electron.vx, electron.vy

        # Current densities
        Jx = q * vx / Vpix
        Jy = q * vy / Vpix

        # Accumulate spatial maps
        self.current_density_map_x[iy, ix] += Jx
        self.current_density_map_y[iy, ix] += Jy
        self.electron_count_map[iy, ix] += 1

        # Continuous injection: offset the timeframe index
        frame = (electron.first_timestep + timestep_number) // self.timesteps_per_timeframe
        if frame < cf.number_of_timeframes:
            self.current_density_profile_x[ix, frame] += Jx
            self.current_density_profile_y[iy, frame] += Jy / cf.number_of_particles
            self.density_profile_y[iy, frame] += 1 / (Vpix*cf.number_of_particles)

    def calculate_current_density_magnitude(self):
        self.current_density_magnitude = np.sqrt(
            self.current_density_map_x**2 + self.current_density_map_y**2
        )

    def calculate_conductivity_from_density_gradient(self):
        e = elementary_charge
        dy = cf.length / cf.number_of_pixels_y
        self.conductivity_from_gradient[:, 0] = (
            np.arange(cf.number_of_timeframes) + 0.5
        ) * self.timesteps_per_timeframe * cf.timestep * 1e9  # t [ns]

        for t in range(cf.number_of_timeframes):
            dt_frame = self.timesteps_per_timeframe * cf.timestep
            J_avg = np.mean(self.current_density_profile_y[:, t]) / dt_frame
            y_coords = np.arange(cf.number_of_pixels_y) * dy
            n_profile = self.density_profile_y[:, t] / cf.number_of_particles

            if np.all(n_profile == 0):
                continue

            slope, _ = np.polyfit(y_coords, n_profile, 1)
            grad_n = slope

            if grad_n != 0:
                sigma = -J_avg / (e * grad_n)
                self.conductivity_from_gradient[t, 1] = sigma

        stable = self.conductivity_from_gradient[cf.number_of_stabilization_timeframes:, 1]
        self.avg_sigma_from_gradient = np.mean(stable)
        self.std_sigma_from_gradient = np.std(stable)

    def write_into_files(self):
        coords_x = np.arange(cf.number_of_pixels_x) * 1e6 * cf.width / cf.number_of_pixels_x
        coords_y = np.arange(cf.number_of_pixels_y) * 1e6 * cf.length / cf.number_of_pixels_y

        data_Jx = np.vstack((coords_x, self.current_density_profile_x.T)).T
        data_Jy = np.vstack((coords_y, self.current_density_profile_y.T)).T
        data_ny = np.vstack((coords_y, self.density_profile_y.T)).T
        
        # Headers for profiles
        hdr_Jx = 'x (m-6),' + ','.join(f'Jx_frame{t+1}' for t in range(cf.number_of_timeframes))
        hdr_Jy = 'y (m-6),' + ','.join(f'Jy_frame{t+1}' for t in range(cf.number_of_timeframes))
        hdr_ny = 'y (m-6),' + ','.join(f'n_frame{t+1}' for t in range(cf.number_of_timeframes))
        hdr_map = ','.join(f'{x:.3f}' for x in coords_x)

        np.savetxt("Data/Current density profile x.csv", data_Jx, fmt='%1.3e', delimiter=",", header=hdr_Jx)
        np.savetxt("Data/Current density profile y.csv", data_Jy, fmt='%1.3e', delimiter=",", header=hdr_Jy)
        np.savetxt("Data/Density profile y.csv", data_ny, fmt='%1.3e', delimiter=",", header=hdr_ny)
        np.savetxt("Data/Current density x.csv", self.current_density_map_x, fmt='%1.2e', delimiter=",", header=hdr_map)
        np.savetxt("Data/Current density y.csv", self.current_density_map_y, fmt='%1.2e', delimiter=",", header=hdr_map)
        np.savetxt("Data/Current density magnitude.csv", self.current_density_magnitude, fmt='%1.2e', delimiter=",", header=hdr_map)
        np.savetxt("Data/Electron count map.csv", self.electron_count_map, fmt='%1.0f', delimiter=",", header=hdr_map)
        np.savetxt("Data/Pixel volumes.csv", self.vol_pixel_ratio, fmt='%1i', delimiter=",", header=hdr_map)

        data_sigma_grad = np.vstack((
            self.conductivity_from_gradient.T,
            np.full(cf.number_of_timeframes, self.avg_sigma_from_gradient),
            np.full(cf.number_of_timeframes, self.std_sigma_from_gradient)
        )).T
        np.savetxt("Data/Conductivity from gradient.csv", data_sigma_grad, fmt='%1.3e', delimiter=",",
                   header="t (ns), sigma (S/m), avg sigma (S/m), std sigma (S/m)", encoding='utf-8')

    def dump_data(self):
        return {
            'current_density_map_x': self.current_density_map_x,
            'current_density_map_y': self.current_density_map_y,
            'current_density_magnitude': self.current_density_magnitude,
            'current_density_profile_x': self.current_density_profile_x,
            'current_density_profile_y': self.current_density_profile_y,
            'density_profile_y': self.density_profile_y,
            'conductivity_from_gradient': self.conductivity_from_gradient,
            'avg_sigma_from_gradient': self.avg_sigma_from_gradient,
            'std_sigma_from_gradient': self.std_sigma_from_gradient,
        }
