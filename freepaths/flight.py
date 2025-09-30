"""Module that provides particle flight characteristics"""

from math import cos
import numpy as np


class Path:
    """Particle path coordinates"""

    def __init__(self, x_init, y_init, z_init):
        """Initialize a path"""
        self.x = np.array([x_init])
        self.y = np.array([y_init])
        self.z = np.array([z_init])

    def add_point(self, x_new, y_new, z_new):
        """Add a point to the particle path"""
        self.x = np.append(self.x, x_new)
        self.y = np.append(self.y, y_new)
        self.z = np.append(self.z, z_new)

    @property
    def number_of_path_points(self):
        """Calculate number of points in the path"""
        return len(self.x)


class Flight:
    """Particle flight characteristics"""

    def __init__(self, particle):
        """Initialize a particle
        flight"""
        self.particle = particle
        self.initial_frequency = self.particle.f
        self.initial_theta = self.particle.theta
        self.path = Path(self.particle.x, self.particle.y, self.particle.z)
        self.exit_theta = 0.0
        self.free_path = 0.0
        self.travel_time = 0.0
        self.time_since_previous_scattering = 0.0
        self.free_paths = []
        self.hole_diff_scattering_angles = []
        self.hole_spec_scattering_angles = []
        self.thermal_conductivity = 0.0
        self.interfaces_angles = []
        self.interfaces_transmission_factor = []
        self.interfaces_wavelength = []
        self.interfaces_frequency = []
        self.interfaces_mode = []

    @property
    def mean_free_path(self):
        """Mean value of all free flights"""
        try:
            mfp = sum(self.free_paths)/len(self.free_paths)
        except:
            mfp = 0
        return mfp

    def add_point_to_path(self):
        """Add a scattering point to the path"""
        self.path.add_point(self.particle.x, self.particle.y, self.particle.z)

    def save_free_paths(self):
        """Save current free path to the list of free paths"""
        self.free_paths.append(self.free_path)

    def save_hole_diff_scattering_angle(self, angle):
        """Save angle of diffuse scattering from the hole"""
        self.hole_diff_scattering_angles.append(angle)

    def save_hole_spec_scattering_angle(self, angle):
        """Save angle of specular scattering from the hole"""
        self.hole_spec_scattering_angles.append(angle)

    def save_interfaces_angles(self, angle):
        """Save a transmission event angle"""
        self.interfaces_angles.append(angle)

    def save_interfaces_transmission_factor(self, alpha_total):
        """Save a transmission event"""
        self.interfaces_transmission_factor.append(alpha_total)

    def save_interfaces_wavelength(self):
        """Save a wavelength event"""
        self.interfaces_wavelength.append(self.particle.wavelength)

    def save_interfaces_frequency(self):
        """Save a frequency event"""
        self.interfaces_frequency.append(self.particle.f)

    def save_interfaces_mode (self):
        """Save mode"""
        self.interfaces_mode.append(self.particle.branch_number)

    def restart(self):
        """Restart the flight after a scattering event"""
        self.time_since_previous_scattering = 0.0
        self.free_path = 0.0

    def finish(self, step, timestep):
        """Finish the flight and record final state"""
        self.exit_theta = self.particle.theta
        self.travel_time = step * timestep

    def add_step(self, timestep):
        """Increase parameters of the flight by length of one step"""
        step_length = self.particle.speed * timestep
        self.free_path += step_length
        self.time_since_previous_scattering += timestep

    def reset_travel_time(self):
        """Reset the travel time of a particle"""
        self.travel_time = 0.0
