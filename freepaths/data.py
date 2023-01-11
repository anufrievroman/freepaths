"""Module that controles recording of various data"""

import numpy as np

from parameters import *

class PathData:
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
        np.savetxt(filename, data, fmt='%1.2e', delimiter=",", header="X (μm), Y (μm), Z (μm)")


class GeneralData:
    """General statistics of various phonon properties"""

    def __init__(self):
        """Initialize arrays for writing various properties"""
        self.initial_angles = []
        self.exit_angles = []
        self.free_paths = []
        self.free_paths_along_y = []
        self.frequencies = []
        self.detected_frequencies = []
        self.group_velocities = []
        self.travel_times = []

    def save_phonon_data(self, ph):
        """Add information about the phonon to the dataset"""
        self.frequencies.append(ph.f)
        self.group_velocities.append(ph.speed)

    def save_flight_data(self, flight):
        """Add information about the phonon flight to the dataset"""
        self.initial_angles.append(flight.initial_theta)
        self.exit_angles.append(flight.exit_theta)
        self.free_paths.extend(flight.free_paths)
        self.free_paths_along_y.extend(flight.free_paths_along_y)
        self.travel_times.append(flight.travel_time)
        self.detected_frequencies.append(flight.detected_frequency)

    def write_into_files(self):
        """Write all the data into files"""
        np.savetxt("Data/All free paths.csv", self.free_paths, fmt='%1.3e', delimiter=",", header="L [m]")
        np.savetxt("Data/All free paths in plane.csv", self.free_paths_along_y, fmt='%1.3e', delimiter=",", header="Ly [m]")
        np.savetxt("Data/All initial frequencies.csv", self.frequencies, fmt='%1.3e', delimiter=",", header="f [Hz]")
        np.savetxt("Data/All detected frequencies.csv", self.detected_frequencies, fmt='%1.3e', delimiter=",", header="f [Hz]")
        np.savetxt("Data/All exit angles.csv", self.exit_angles, fmt='%1.3e', delimiter=",", header="Angle [rad]")
        np.savetxt("Data/All initial angles.csv", self.initial_angles, fmt='%1.3e', delimiter=",", header="Angle [rad]")
        np.savetxt("Data/All group velocities.csv", self.group_velocities, fmt='%1.3e', delimiter=",", header="Vg [rad]")
        np.savetxt("Data/All travel times.csv", self.travel_times, fmt='%1.3e', delimiter=",", header="Travel time [s]")


class ScatteringData:
    """Statistics of phonon scattering events"""

    def __init__(self):
        """Initialize arrays according to the number of segments"""
        self.wall_diffuse = np.zeros(NUMBER_OF_LENGTH_SEGMENTS)
        self.wall_specular = np.zeros(NUMBER_OF_LENGTH_SEGMENTS)
        self.top_diffuse = np.zeros(NUMBER_OF_LENGTH_SEGMENTS)
        self.top_specular = np.zeros(NUMBER_OF_LENGTH_SEGMENTS)
        self.hole_diffuse = np.zeros(NUMBER_OF_LENGTH_SEGMENTS)
        self.hole_specular = np.zeros(NUMBER_OF_LENGTH_SEGMENTS)
        self.pillar_diffuse = np.zeros(NUMBER_OF_LENGTH_SEGMENTS)
        self.pillar_specular = np.zeros(NUMBER_OF_LENGTH_SEGMENTS)
        self.hot_side = np.zeros(NUMBER_OF_LENGTH_SEGMENTS)
        self.internal = np.zeros(NUMBER_OF_LENGTH_SEGMENTS)
        self.total = np.zeros(NUMBER_OF_LENGTH_SEGMENTS)

    def save_scattering_events(self, y, scattering_types):
        """Analyze types of scattering at the current timestep and add it to the statistics"""

        # Calculate in which length segment (starting from zero) we are:
        segment = int(y // (LENGTH / NUMBER_OF_LENGTH_SEGMENTS))
        self.total[segment] += 1

        # Scattering on side walls:
        self.wall_diffuse[segment]  += 1 if scattering_types.walls == 'diffuse' else 0
        self.wall_specular[segment] += 1 if scattering_types.walls == 'specular' else 0

        # Scattering on top and bottom:
        self.top_diffuse[segment]  += 1 if scattering_types.top_bottom == 'diffuse' else 0
        self.top_specular[segment] += 1 if scattering_types.top_bottom == 'specular' else 0

        # Scattering on holes:
        self.hole_diffuse[segment]  += 1 if scattering_types.holes == 'diffuse' else 0
        self.hole_specular[segment] += 1 if scattering_types.holes == 'specular' else 0

        # Scattering on pillars:
        self.pillar_diffuse[segment]  += 1 if scattering_types.pillars == 'diffuse' else 0
        self.pillar_specular[segment] += 1 if scattering_types.pillars == 'specular' else 0

        # Internal scattering and rethermalization on hot side:
        self.hot_side[segment] += 1 if scattering_types.hot_side == 'diffuse' else 0
        self.internal[segment] += 1 if scattering_types.internal == 'diffuse' else 0

    def write_into_files(self):
        """Write data into a file"""
        filename = "Data/Scattering events statistics.csv"
        data = np.vstack((self.wall_diffuse, self.wall_specular, self.top_diffuse, self.top_specular, self.hole_diffuse,
                        self.hole_specular, self.hot_side, self.internal, self.pillar_diffuse, self.pillar_specular)).T
        header1 = "Sidewalls diffuse, Sidewalls specular, Top & bottom diffuse, Top & bottom specular, "
        header2 = "Holes diffuse, Holes specular, Hot side, Internal, Pillars diffuse, Pillars specular"
        header = header1 + header2
        np.savetxt(filename, data, fmt='%1.3e', delimiter=",", header=header)


class SegmentData:
    """Statistics of events happening in different segments"""

    def __init__(self):
        """Initialize arrays according to the number of segments"""
        self.time_spent = np.zeros(NUMBER_OF_LENGTH_SEGMENTS)

    @property
    def segment_coordinates(self):
        """Calculate coordinates of the centers of each segment"""
        segment_length = LENGTH / NUMBER_OF_LENGTH_SEGMENTS
        segments = [(segment_length/2 + i*segment_length) * 1e6 for i in range(NUMBER_OF_LENGTH_SEGMENTS)]
        return segments

    def record_time_in_segment(self, coordinate):
        """Record how long phonon stays in different segments"""
        for segment_number in range(NUMBER_OF_LENGTH_SEGMENTS):
            segment_beginning = segment_number * (LENGTH / NUMBER_OF_LENGTH_SEGMENTS)
            segment_end = (segment_number + 1)*(LENGTH / NUMBER_OF_LENGTH_SEGMENTS)
            if segment_beginning <= coordinate < segment_end:
                self.time_spent[segment_number] += TIMESTEP * 1e6

    def write_into_files(self):
        """Write data into files"""
        filename = "Data/Time spent in segments.csv"
        data = np.vstack((self.segment_coordinates, self.time_spent)).T
        np.savetxt(filename, data, fmt='%1.3e', delimiter=",", header="Y [um], Time [ns]")