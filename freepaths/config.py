"""Module that reads the user input file, privides default values, and converts the variables into enums"""

import sys
import argparse

from freepaths.options import Materials, Distributions, Positions

# Import a default input file:
from freepaths.default_config import *


# Parse user arguments:
parser = argparse.ArgumentParser(prog = 'FreePATHS', description = 'Monte Carlo simulator',
                epilog = 'For more information, visit: https://github.com/anufrievroman/freepaths')
parser.add_argument('input_file', nargs='?', default=None, help='The input file')
parser.add_argument("-s", "--sampling", help="Run in MFP sampling mode", action="store_true")
args = parser.parse_args()


# If a file is provided, overwrite the default values:
if args.input_file:
    exec(open(args.input_file).read(), globals())
else:
    print("You didn't provide any input file, so we run a demo simulation:\n")


class Config:
    """Class that contains all the settings for the simulation"""

    def __init__(self):
        """Initiate all the parameters from global variables"""

        # General parameters:
        self.output_folder_name = OUTPUT_FOLDER_NAME
        self.number_of_phonons = NUMBER_OF_PHONONS
        self.number_of_timesteps = NUMBER_OF_TIMESTEPS
        self.number_of_nodes = NUMBER_OF_NODES
        self.timestep = TIMESTEP
        self.temp = T
        self.plots_in_terminal = PLOTS_IN_TERMINAL
        self.output_scattering_map = OUTPUT_SCATTERING_MAP
        self.output_raw_thermal_map = OUTPUT_RAW_THERMAL_MAP
        self.output_trajectories_of_first = OUTPUT_TRAJECTORIES_OF_FIRST
        self.number_of_length_segments = NUMBER_OF_LENGTH_SEGMENTS
        self.hot_side_angle_distribution = HOT_SIDE_ANGLE_DISTRIBUTION

        # Map & profiles parameters:
        self.number_of_pixels_x = NUMBER_OF_PIXELS_X
        self.number_of_pixels_y = NUMBER_OF_PIXELS_Y
        self.number_of_timeframes = NUMBER_OF_TIMEFRAMES

        # Material parameters:
        self.media = MEDIA
        self.specific_heat_capacity = SPECIFIC_HEAT_CAPACITY

        # Internal scattering:
        self.include_internal_scattering = INCLUDE_INTERNAL_SCATTERING
        self.use_gray_approximation_mfp = USE_GRAY_APPROXIMATION_MFP
        self.gray_approximation_mfp = GRAY_APPROXIMATION_MFP

        # System dimensions:
        self.thickness = THICKNESS
        self.width = WIDTH
        self.length = LENGTH

        # Hot and cold sides:
        self.frequency_detector_size = FREQUENCY_DETECTOR_SIZE
        self.cold_side_position = COLD_SIDE_POSITION
        self.hot_size_x = HOT_SIZE_X
        self.hot_size_width = HOT_SIZE_WIDTH

        # Roughness:
        self.side_wall_roughness = SIDE_WALL_ROUGHNESS
        self.hole_roughness = HOLE_ROUGHNESS
        self.pillar_roughness = PILLAR_ROUGHNESS
        self.top_roughness = TOP_ROUGHNESS
        self.bottom_roughness = BOTTOM_ROUGHNESS
        self.pillar_top_roughness = PILLAR_TOP_ROUGHNESS

        # Hole array parameters:
        self.include_holes = INCLUDE_HOLES
        self.circular_hole_diameter = CIRCULAR_HOLE_DIAMETER
        self.rectangular_hole_side_x = RECTANGULAR_HOLE_SIDE_X
        self.rectangular_hole_side_y = RECTANGULAR_HOLE_SIDE_Y
        self.period_x = PERIOD_X
        self.period_y = PERIOD_Y

        # Lattice of holes:
        self.hole_coordinates = HOLE_COORDINATES
        self.hole_shapes = HOLE_SHAPES

        # Pillar array parameters [m]
        self.include_pillars = INCLUDE_PILLARS
        self.pillar_coordinates = PILLAR_COORDINATES
        self.pillar_height = PILLAR_HEIGHT
        self.pillar_wall_angle = PILLAR_WALL_ANGLE


    def convert_to_enums(self):
        """Convert some user generated parameters into enums"""
        if self.hot_side_angle_distribution in ["random", "lambert", "directional"]:
            self.hot_side_angle_distribution = Distributions[self.hot_side_angle_distribution.upper()]
        else:
            print(f"ERROR: Parameter {self.hot_side_angle_distribution} is not set correctly.\n")

        if self.media in ["Si", "SiC", "Diamond", "AlN"]:
            self.media = Materials[self.media]
        else:
            print(f"ERROR: Material {self.media} is not in the database.\n")

        if self.cold_side_position in ["top", "top_and_right", "top_and_bottom", "right"]:
            self.cold_side_position = Positions[self.cold_side_position.upper()]
        else:
            print(f"ERROR: Parameter {self.cold_side_position} is not set correctly.\n")


    def check_validity(self):
        """Check if some of the parameteres are valid"""
        if self.number_of_phonons < self.output_trajectories_of_first:
            print("ERROR: Parameter OUTPUT_TRAJECTORIES_OF_FIRST exeeds NUMBER_OF_PHONONS!\n")


cf = Config()
cf.convert_to_enums()
cf.check_validity()