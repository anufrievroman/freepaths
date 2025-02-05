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
    """Paths of phonons in space"""

    def __init__(self):
        """Initialize arrays to save phonon paths"""
        self.phonon_paths = []

    def save_phonon_path(self, flight):
        """Save the path to list of all paths"""
        self.phonon_paths.append(flight.path)

    @property
    def length_of_longest_path(self):
        """Calculate the number of points in the longest path"""
        return max([path.number_of_path_points for path in self.phonon_paths])

    def write_into_files(self):
        """Write all the path coordinates into a file"""
        filename = "Data/Phonon paths.csv"
        data = np.zeros((self.length_of_longest_path, 3*len(self.phonon_paths)))
        for index, path in enumerate(self.phonon_paths):
            for point_n, (x, y, z) in enumerate(zip(path.x, path.y, path.z)):
                data[point_n, 0 + index*3] = x*1e6
                data[point_n, 1 + index*3] = y*1e6
                data[point_n, 2 + index*3] = z*1e6
        np.savetxt(filename, data, fmt='%2.4f', delimiter=",", header="X (μm), Y (μm), Z (μm)", encoding='utf-8')

    def dump_data(self):
        """Return data of a process in the form of a dictionary to be attached to the global data"""
        return {'phonon_paths': self.phonon_paths}


class GeneralData(Data):
    """General statistics of various phonon properties"""

    def __init__(self):
        """Initialize arrays for writing various properties"""
        self.initial_angles = []
        self.exit_angles = []
        self.hole_diff_scattering_angles = []
        self.hole_spec_scattering_angles = []
        self.free_paths = []
        self.free_paths_along_y = []
        self.frequencies = []
        self.group_velocities = []
        self.travel_times = []
        self.mean_free_paths = []
        self.thermal_conductivity = []

    def save_phonon_data(self, ph):
        """Add information about the phonon to the dataset"""
        self.frequencies.append(ph.f)
        self.group_velocities.append(ph.speed)

    def save_flight_data(self, flight):
        """Add information about the phonon flight to the dataset"""
        self.initial_angles.append(flight.initial_theta)
        self.exit_angles.append(flight.exit_theta)
        self.free_paths.extend(flight.free_paths)
        self.hole_diff_scattering_angles.extend(flight.hole_diff_scattering_angles)
        self.hole_spec_scattering_angles.extend(flight.hole_spec_scattering_angles)
        self.free_paths_along_y.extend(flight.free_paths_along_y)
        self.travel_times.append(flight.travel_time)
        self.mean_free_paths.append(flight.mean_free_path)
        self.thermal_conductivity.append(flight.thermal_conductivity)

    def write_into_files(self):
        """Write all the data into files"""
        np.savetxt("Data/All free paths.csv", self.free_paths, fmt='%2.4e', header="L [m]", encoding='utf-8')
        np.savetxt("Data/All free paths in plane.csv", self.free_paths_along_y, fmt='%2.4e', header="Ly [m]", encoding='utf-8')
        np.savetxt("Data/All initial frequencies.csv", self.frequencies, fmt='%2.4e', header="f [Hz]", encoding='utf-8')
        np.savetxt("Data/All exit angles.csv", self.exit_angles, fmt='%.4f', header="Angle [rad]", encoding='utf-8')
        np.savetxt("Data/All hole diffuse scattering angles.csv", self.hole_diff_scattering_angles, fmt='%.4f', header="Angle [rad]", encoding='utf-8')
        np.savetxt("Data/All hole specular scattering angles.csv", self.hole_spec_scattering_angles, fmt='%.4f', header="Angle [rad]", encoding='utf-8')
        np.savetxt("Data/All initial angles.csv", self.initial_angles, fmt='%.4f', header="Angle [rad]", encoding='utf-8')
        np.savetxt("Data/All group velocities.csv", self.group_velocities, fmt='%.4f', header="Vg [m//s]", encoding='utf-8')
        np.savetxt("Data/All travel times.csv", self.travel_times, fmt='%2.4e', header="Travel time [s]", encoding='utf-8')
        np.savetxt("Data/All mean free paths.csv", self.mean_free_paths, fmt='%2.4e', header="MFPs [m]", encoding='utf-8')
        np.savetxt("Data/All thermal conductivities.csv", self.thermal_conductivity, fmt='%2.4e', header="K [W/mK]", encoding='utf-8')

    def dump_data(self):
        """Return data of a process in the form of a dictionary to be attached to the global data"""
        return {
            'initial_angles': self.initial_angles,
            'exit_angles': self.exit_angles,
            'hole_diff_scattering_angles': self.hole_diff_scattering_angles,
            'hole_spec_scattering_angles': self.hole_spec_scattering_angles,
            'free_paths': self.free_paths,
            'free_paths_along_y': self.free_paths_along_y,
            'frequencies': self.frequencies,
            'group_velocities': self.group_velocities,
            'travel_times': self.travel_times,
            'mean_free_paths': self.mean_free_paths,
            'thermal_conductivity': self.thermal_conductivity,
        }


class ScatteringData(Data):
    """Statistics of phonon scattering events"""

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
        except:
            pass

    def write_into_files(self):
        """Write data into a file"""
        filename = "Data/Scattering events statistics.csv"
        data = np.vstack((self.wall_diffuse, self.wall_specular, self.top_diffuse, self.top_specular, self.hole_diffuse,
                        self.hole_specular, self.hot_side, self.internal, self.pillar_diffuse, self.pillar_specular)).T
        header1 = "Sidewalls diffuse, Sidewalls specular, Top & bottom diffuse, Top & bottom specular, "
        header2 = "Holes diffuse, Holes specular, Hot side, Internal, Pillars diffuse, Pillars specular"
        header = header1 + header2
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
            'total': self.total
        }


class TriangleScatteringData(Data):
    """Statistics of phonon scattering events on triangular holes"""

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
        """Record how long phonon stays in different segments"""
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

