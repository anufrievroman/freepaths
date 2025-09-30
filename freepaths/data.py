"""Module that controls recording various data"""

import numpy as np

from freepaths.config import cf
from freepaths.scattering_types import Scattering

class Data:
    """Parent data class with functions common to all classes below"""

    def read_data(self, data_dict):
        """Read the data from the finished worker and add new data to already existing"""
        for key, value in data_dict.items():
            setattr(self, key, getattr(self, key) + value)


class PathData(Data):
    """Paths of particles in space"""

    def __init__(self):
        """Initialize arrays to save particle paths"""
        self.particle_paths = []

    def save_particle_path(self, flight):
        """Save the path to list of all paths"""
        if cf.low_memory_usage:
            return
        self.particle_paths.append(flight.path)

    @property
    def length_of_longest_path(self):
        """Calculate the number of points in the longest path"""
        return max([path.number_of_path_points for path in self.particle_paths])

    def write_into_files(self):
        """Write all the path coordinates into a file"""
        if cf.low_memory_usage:
            return
        filename = "Data/Particle paths.csv"
        data = np.zeros((self.length_of_longest_path, 3*len(self.particle_paths)))
        for index, path in enumerate(self.particle_paths):
            for point_n, (x, y, z) in enumerate(zip(path.x, path.y, path.z)):
                data[point_n, 0 + index*3] = x*1e6
                data[point_n, 1 + index*3] = y*1e6
                data[point_n, 2 + index*3] = z*1e6
        np.savetxt(filename, data, fmt='%2.4f', delimiter=",", header="X (μm), Y (μm), Z (μm)", encoding='utf-8')

    def dump_data(self):
        """Return data of a process in the form of a dictionary to be attached to the global data"""
        return {'particle_paths': self.particle_paths}


class GeneralData(Data):
    """General statistics of various particle properties"""

    def __init__(self):
        """Initialize arrays for writing various properties"""
        self.initial_angles = []
        self.exit_angles = []
        self.hole_diff_scattering_angles = []
        self.hole_spec_scattering_angles = []
        self.free_paths = []
        self.frequencies = []
        self.initial_energies = []
        self.group_velocities = []
        self.travel_times = []
        self.mean_free_paths = []
        self.thermal_conductivity = []
        self.interfaces_angles = []
        self.interfaces_transmission_factor = []
        self.interfaces_wavelength = []
        self.interfaces_frequency = []
        self.interfaces_mode = []

    def save_particle_data(self, pt):
        """Add information about the particle to the dataset"""
        self.frequencies.append(pt.f)
        self.group_velocities.append(pt.speed)
        self.initial_energies.append(pt.energy)

    def save_flight_data(self, flight):
        """Add information about the particle flight to the dataset"""
        if not cf.low_memory_usage:
            self.free_paths.extend(flight.free_paths)
        self.initial_angles.append(flight.initial_theta)
        self.exit_angles.append(flight.exit_theta)
        self.hole_diff_scattering_angles.extend(flight.hole_diff_scattering_angles)
        self.hole_spec_scattering_angles.extend(flight.hole_spec_scattering_angles)
        self.travel_times.append(flight.travel_time)
        self.mean_free_paths.append(flight.mean_free_path)
        self.thermal_conductivity.append(flight.thermal_conductivity)
        self.interfaces_angles.extend(flight.interfaces_angles)
        self.interfaces_transmission_factor.extend(flight.interfaces_transmission_factor)
        self.interfaces_wavelength.extend(flight.interfaces_wavelength)
        self.interfaces_frequency.extend(flight.interfaces_frequency)
        self.interfaces_mode.extend(flight.interfaces_mode)



    def write_into_files(self):
        import os
        if not os.path.exists("Data"):
            os.makedirs("Data")
        """Write all the data into files"""
        if not cf.low_memory_usage:
            np.savetxt("Data/All free paths.csv", self.free_paths, fmt='%2.4e', header="L [m]", encoding='utf-8')
        np.savetxt("Data/All initial frequencies.csv", self.frequencies, fmt='%2.4e', header="f [Hz]", encoding='utf-8')
        np.savetxt("Data/All exit angles.csv", self.exit_angles, fmt='%.4f', header="Angle [rad]", encoding='utf-8')
        np.savetxt("Data/All hole diffuse scattering angles.csv", self.hole_diff_scattering_angles, fmt='%.4f', header="Angle [rad]", encoding='utf-8')
        np.savetxt("Data/All hole specular scattering angles.csv", self.hole_spec_scattering_angles, fmt='%.4f', header="Angle [rad]", encoding='utf-8')
        np.savetxt("Data/All initial angles.csv", self.initial_angles, fmt='%.4f', header="Angle [rad]", encoding='utf-8')
        np.savetxt("Data/All group velocities.csv", self.group_velocities, fmt='%.4f', header="Vg [m//s]", encoding='utf-8')
        np.savetxt("Data/All travel times.csv", self.travel_times, fmt='%2.4e', header="Travel time [s]", encoding='utf-8')
        np.savetxt("Data/All mean free paths.csv", self.mean_free_paths, fmt='%2.4e', header="MFPs [m]", encoding='utf-8')
        np.savetxt("Data/All thermal conductivities.csv", self.thermal_conductivity, fmt='%2.4e', header="K [W/mK]", encoding='utf-8')
        np.savetxt("Data/All initial energies.csv", self.initial_energies, fmt='%2.4e', header="Energy [J]", encoding='utf-8')
        np.savetxt("Data/All interfaces angles.csv", self.interfaces_angles, fmt='%.4f', header="Angle [rad]", encoding='utf-8')
        np.savetxt("Data/All interfaces transmission factor.csv", self.interfaces_transmission_factor, fmt='%2.4e', header="Transmission factor [%]", encoding='utf-8')
        np.savetxt("Data/All interfaces wavelength.csv", self.interfaces_wavelength, fmt='%.18f', header="Interfaces wavelength [m]", encoding='utf-8')
        np.savetxt("Data/All interfaces frequency.csv", self.interfaces_frequency, fmt='%.18f', header="Interfaces frequency [THz]", encoding='utf-8')
        np.savetxt("Data/All interfaces mode.csv", self.interfaces_mode, fmt='%d', header="Interfaces mode", encoding='utf-8')

    def dump_data(self):
        """Return data of a process in the form of a dictionary to be attached to the global data"""
        return {
            'initial_angles': self.initial_angles,
            'exit_angles': self.exit_angles,
            'hole_diff_scattering_angles': self.hole_diff_scattering_angles,
            'hole_spec_scattering_angles': self.hole_spec_scattering_angles,
            'free_paths': self.free_paths,
            'frequencies': self.frequencies,
            'group_velocities': self.group_velocities,
            'travel_times': self.travel_times,
            'mean_free_paths': self.mean_free_paths,
            'thermal_conductivity': self.thermal_conductivity,
            'initial_energies': self.initial_energies,
            'interfaces_angles': self.interfaces_angles,
            'interfaces_transmission_factor': self.interfaces_transmission_factor,
            'interfaces_wavelength': self.interfaces_wavelength,
            'interfaces_frequency': self.interfaces_frequency,
            'interfaces_mode': self.interfaces_mode,
        }



class ScatteringData(Data):
    """Statistics of particle scattering events"""

    def __init__(self):
        """Initialize arrays according to the number of segments"""
        self.wall_diffuse = np.zeros(cf.number_of_length_segments+1)
        self.wall_specular = np.zeros(cf.number_of_length_segments+1)
        self.top_diffuse = np.zeros(cf.number_of_length_segments+1)
        self.top_specular = np.zeros(cf.number_of_length_segments+1)
        self.hole_diffuse = np.zeros(cf.number_of_length_segments+1)
        self.hole_specular = np.zeros(cf.number_of_length_segments+1)
        self.pillar_diffuse = np.zeros(cf.number_of_length_segments+1)
        self.pillar_specular = np.zeros(cf.number_of_length_segments+1)
        self.hot_side = np.zeros(cf.number_of_length_segments+1)
        self.internal = np.zeros(cf.number_of_length_segments+1)
        self.interfaces_diffuse = np.zeros(cf.number_of_length_segments+1)
        self.interfaces_specular = np.zeros(cf.number_of_length_segments+1)
        self.total = np.zeros(cf.number_of_length_segments+1)
        self.interfaces_transmission_specular= np.zeros(cf.number_of_length_segments + 1)
        self.interfaces_transmission_diffuse = np.zeros(cf.number_of_length_segments + 1)
        self.interfaces_angles = np.zeros(cf.number_of_length_segments + 1)
        self.interfaces_transmission_factor = np.zeros(cf.number_of_length_segments + 1)
        self.interfaces_wavelength = np.zeros(cf.number_of_length_segments+1)
        self.interfaces_frequency = np.zeros(cf.number_of_length_segments+1)
        self.interfaces_mode = np.zeros(cf.number_of_length_segments+1)



    def save_scattering_events(self, y, scattering_types):
        """Analyze types of scattering at the current timestep and add it to the statistics"""

        try:
            # Calculate in which length segment (starting from zero) we are:
            segment = int(y // (cf.length / cf.number_of_length_segments))
            self.total[segment] += 1

            # Scattering on side walls:
            self.wall_diffuse[segment]  += 1 if scattering_types.walls == Scattering.DIFFUSE else 0
            self.wall_specular[segment] += 1 if scattering_types.walls == Scattering.SPECULAR else 0

            # Scattering on top and bottom:
            self.top_diffuse[segment]  += 1 if scattering_types.top_bottom == Scattering.DIFFUSE else 0
            self.top_specular[segment] += 1 if scattering_types.top_bottom == Scattering.SPECULAR else 0

            # Scattering on holes:
            self.hole_diffuse[segment]  += 1 if scattering_types.holes == Scattering.DIFFUSE else 0
            self.hole_specular[segment] += 1 if scattering_types.holes == Scattering.SPECULAR else 0

            # Scattering on pillars:
            self.pillar_diffuse[segment]  += 1 if scattering_types.pillars == Scattering.DIFFUSE else 0
            self.pillar_specular[segment] += 1 if scattering_types.pillars == Scattering.SPECULAR else 0

            # Scattering on pillars:
            self.interfaces_diffuse[segment]  += 1 if scattering_types.interfaces == Scattering.DIFFUSE else 0
            self.interfaces_specular[segment] += 1 if scattering_types.interfaces == Scattering.SPECULAR else 0

            # Internal scattering and rethermalization on hot side:
            self.hot_side[segment] += 1 if scattering_types.hot_side == Scattering.DIFFUSE else 0
            self.internal[segment] += 1 if scattering_types.internal == Scattering.DIFFUSE else 0

            # Interfaces transmission events:
            self.interfaces_transmission_specular[segment] += 1 if scattering_types.interfaces_transmission== Scattering.SPECULAR else 0
            self.interfaces_transmission_diffuse[segment] += 1 if scattering_types.interfaces_transmission == Scattering.DIFFUSE else 0
        except:
            pass

    def write_into_files(self):
        """Write data into a file"""
        filename = "Data/Scattering events statistics.csv"
        data = np.vstack((self.wall_diffuse, self.wall_specular, self.top_diffuse, self.top_specular, self.hole_diffuse,
                        self.hole_specular, self.hot_side, self.internal, self.pillar_diffuse, self.pillar_specular,
                        self.interfaces_diffuse, self.interfaces_specular,
                        self.interfaces_transmission_diffuse, self.interfaces_transmission_specular,
                        self.interfaces_transmission_factor, self.interfaces_angles,
                        self.interfaces_wavelength, self.interfaces_mode)).T
        header1 = "Sidewalls diffuse, Sidewalls specular, Top & bottom diffuse, Top & bottom specular, "
        header2 = "Holes diffuse, Holes specular, Hot side, Internal, Pillars diffuse, Pillars specular"
        header3 = "Interfaces diffuse, Interfaces specular, Interfaces transmission diffuse, Interfaces transmission specular, " \
        "Interfaces transmission factor, Interfaces angles, Interfaces wavelength, Interfaces frequency, Interfaces mode"
        header = header1 + header2 + header3
        np.savetxt(filename, data, fmt='%1.3e', delimiter=",", header=header, encoding='utf-8')

    def dump_data(self):
        """Return data of a process in the form of a dictionary to be attached to the global data"""
        return {
            'wall_diffuse': self.wall_diffuse,
            'wall_specular': self.wall_specular,
            'top_diffuse': self.top_diffuse,
            'top_specular': self.top_specular,
            'hole_diffuse': self.hole_diffuse,
            'hole_specular': self.hole_specular,
            'pillar_diffuse': self.pillar_diffuse,
            'pillar_specular': self.pillar_specular,
            'hot_side': self.hot_side,
            'internal': self.internal,
            'interfaces_diffuse': self.interfaces_diffuse,
            'interfaces_specular': self.interfaces_specular,
            'interfaces_transmission_specular': self.interfaces_transmission_specular,
            'interfaces_transmission_diffuse': self.interfaces_transmission_diffuse,
            'interfaces_angles': self.interfaces_angles,
            'interfaces_transmission_factor': self.interfaces_transmission_factor,
            'interfaces_wavelength': self.interfaces_wavelength,
            'interfaces_frequency': self.interfaces_frequency,
            'interfaces_mode' : self.interfaces_mode,
            'total': self.total
        }


class TriangleScatteringData(Data):
    """Statistics of particle scattering events on triangular holes"""

    def __init__(self):
        """Initialize arrays according to the number of segments"""
        self.right_wall_diffuse = 0
        self.right_wall_specular = 0
        self.left_wall_diffuse = 0
        self.left_wall_specular = 0
        self.floor_diffuse = 0
        self.floor_specular = 0

    def save_scattering_events(self, y, triangle_scattering_places):
        """Analyze types of scattering at the current timestep and add it to the statistics"""

        try:
            self.right_wall_diffuse  += 1 if triangle_scattering_places.right_wall == Scattering.DIFFUSE else 0
            self.right_wall_specular += 1 if triangle_scattering_places.right_wall == Scattering.SPECULAR else 0
            self.left_wall_diffuse  += 1 if triangle_scattering_places.left_wall == Scattering.DIFFUSE else 0
            self.left_wall_specular += 1 if triangle_scattering_places.left_wall == Scattering.SPECULAR else 0
            self.floor_diffuse += 1 if triangle_scattering_places.floor == Scattering.DIFFUSE else 0
            self.floor_specular += 1 if triangle_scattering_places.floor == Scattering.SPECULAR else 0
        except:
            pass

    def write_into_files(self):
        """Write data into a file"""
        filename = "Data/Triangle scattering statistics.csv"
        data = np.vstack(np.array([self.right_wall_diffuse, self.right_wall_specular,
                self.left_wall_diffuse, self.left_wall_specular,
                self.floor_diffuse, self.floor_specular])).T
        header = "Inclined right diffuse, Inclinded right specular, Inclined left diffuse, Inclinded left specular, Floor diffuse, Floor specular"
        np.savetxt(filename, data, fmt='%1.3e', delimiter=",", header=header, encoding='utf-8')

    def dump_data(self):
        """Return data of a process in the form of a dictionary to be attached to the global data"""
        return {
            'right_wall_diffuse': self.right_wall_diffuse,
            'right_wall_specular': self.right_wall_specular,
            'left_wall_diffuse': self.left_wall_diffuse,
            'left_wall_specular': self.left_wall_specular,
            'floor_diffuse': self.floor_diffuse,
            'floor_specular': self.floor_specular,
        }



class SegmentData(Data):
    """Statistics of events happening in different segments"""

    def __init__(self):
        """Initialize arrays according to the number of segments"""
        self.time_spent = np.zeros(cf.number_of_length_segments)

    @property
    def segment_coordinates(self):
        """Calculate coordinates of the centers of each segment"""
        segment_length = cf.length * 1e6 / cf.number_of_length_segments
        segments = [(segment_length/2 + i*segment_length) for i in range(cf.number_of_length_segments)]
        return segments

    def record_time_in_segment(self, coordinate):
        """Record how long particle stays in different segments"""
        for segment_number in range(cf.number_of_length_segments):
            segment_beginning = segment_number * (cf.length / cf.number_of_length_segments)
            segment_end = (segment_number + 1)*(cf.length / cf.number_of_length_segments)
            if segment_beginning <= coordinate < segment_end:
                self.time_spent[segment_number] += cf.timestep * 1e6

    def write_into_files(self):
        """Write data into files"""
        filename = "Data/Time spent in segments.csv"
        data = np.vstack((self.segment_coordinates, self.time_spent)).T
        np.savetxt(filename, data, fmt='%1.3e', delimiter=",", header="Y [um], Time [ns]", encoding='utf-8')

    def dump_data(self):
        """Return data of a process in the form of a dictionary to be attached to the global data"""
        return {'time_spent': self.time_spent}

