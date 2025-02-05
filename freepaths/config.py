"""
Module that reads the user input file, provides default values,
checks the validity of the parameters and converts the variables into enums
"""

import sys
import argparse
import logging
from colorama import Fore, Style

from freepaths.sources import Distributions
from freepaths.scatterers import *

# Import a default input file:
from freepaths.default_config import *

# Start logging:
logging.basicConfig(level=logging.WARNING, format=f"{Fore.RED}%(levelname)s:{Style.RESET_ALL} %(message)s",
                    handlers=[logging.StreamHandler(), ])

# Parse user arguments:
WEBSITE = 'https://anufrievroman.gitbook.io/freepaths'
parser = argparse.ArgumentParser(prog='FreePATHS', description='Monte Carlo simulator',
                                 epilog=f'For more information, visit: {WEBSITE}')
parser.add_argument('input_file', nargs='?', default=None, help='The input file')
parser.add_argument("-s", "--sampling", help="Run in MFP sampling mode", action="store_true")
args = parser.parse_args()


# If a file is provided, overwrite the default values:
if args.input_file:
    try:
        exec(open(args.input_file, encoding='utf-8').read(), globals())
    except FileNotFoundError:
        logging.error("Input file does not exist. Check if you are in the right folder and the file name is correct.")
        sys.exit()
else:
    logging.warning("You provided no input file, so we will run a demo simulation:")


class Config:
    """Class that contains all the settings for the simulation"""

    def __init__(self):
        """Initiate all the parameters from global variables"""

        # General parameters:
        self.output_folder_name = str(OUTPUT_FOLDER_NAME)
        self.number_of_phonons = NUMBER_OF_PHONONS
        self.number_of_nodes = NUMBER_OF_NODES
        self.temp = T
        self.output_scattering_map = OUTPUT_SCATTERING_MAP
        self.output_trajectories_of_first = OUTPUT_TRAJECTORIES_OF_FIRST
        self.output_structure_color = OUTPUT_STRUCTURE_COLOR
        self.number_of_length_segments = NUMBER_OF_LENGTH_SEGMENTS

        # Time parameters:
        self.timestep = TIMESTEP
        self.number_of_timesteps = NUMBER_OF_TIMESTEPS
        self.number_of_timeframes = NUMBER_OF_TIMEFRAMES
        self.number_of_stabilization_timeframes = NUMBER_OF_STABILIZATION_TIMEFRAMES

        # Animation:
        self.output_path_animation = OUTPUT_PATH_ANIMATION
        self.output_animation_fps = OUTPUT_ANIMATION_FPS

        # Map & profiles parameters:
        self.number_of_pixels_x = NUMBER_OF_PIXELS_X
        self.number_of_pixels_y = NUMBER_OF_PIXELS_Y
        self.number_of_virtual_timesteps = NUMBER_OF_VIRTUAL_TIMESTEPS
        self.ignore_faulty_phonons = IGNORE_FAULTY_PHONONS

        # Material parameters:
        self.media = MEDIA

        # Internal scattering:
        self.include_internal_scattering = INCLUDE_INTERNAL_SCATTERING
        self.use_gray_approximation_mfp = USE_GRAY_APPROXIMATION_MFP
        self.gray_approximation_mfp = GRAY_APPROXIMATION_MFP

        # System dimensions:
        self.thickness = THICKNESS
        self.width = WIDTH
        self.length = LENGTH
        self.is_two_dimensional_material = IS_TWO_DIMENSIONAL_MATERIAL
        self.include_right_sidewall = INCLUDE_RIGHT_SIDEWALL
        self.include_left_sidewall = INCLUDE_LEFT_SIDEWALL
        self.include_top_sidewall = INCLUDE_TOP_SIDEWALL
        self.include_bottom_sidewall = INCLUDE_BOTTOM_SIDEWALL

        # Hot side positions:
        self.hot_side_position_top = HOT_SIDE_POSITION_TOP
        self.hot_side_position_bottom = HOT_SIDE_POSITION_BOTTOM
        self.hot_side_position_right = HOT_SIDE_POSITION_RIGHT
        self.hot_side_position_left = HOT_SIDE_POSITION_LEFT

        # Sources:
        self.phonon_sources = PHONON_SOURCES

        # Cold side positions:
        self.cold_side_position_top = COLD_SIDE_POSITION_TOP
        self.cold_side_position_bottom = COLD_SIDE_POSITION_BOTTOM
        self.cold_side_position_right = COLD_SIDE_POSITION_RIGHT
        self.cold_side_position_left = COLD_SIDE_POSITION_LEFT

        # Roughness:
        self.side_wall_roughness = SIDE_WALL_ROUGHNESS
        self.hole_roughness = HOLE_ROUGHNESS
        self.pillar_roughness = PILLAR_ROUGHNESS
        self.top_roughness = TOP_ROUGHNESS
        self.bottom_roughness = BOTTOM_ROUGHNESS
        self.pillar_top_roughness = PILLAR_TOP_ROUGHNESS
        self.interface_roughness = INTERFACE_ROUGHNESS

        # Hole array parameters:
        self.holes = HOLES
        self.pillars = PILLARS
        self.interfaces = INTERFACES

        # Multiprocessing:
        self.num_workers = NUMBER_OF_PROCESSES

    def convert_to_enums(self):
        """Convert some user generated parameters into enums"""

        # Distributions:
        valid_distributions =[member.name.lower() for member in Distributions]
        for source in self.phonon_sources:
            if source.angle_distribution in valid_distributions:
                source.angle_distribution = Distributions[source.angle_distribution.upper()]
            else:
                logging.error("Parameter angle_distribution of a source is not set correctly.")
                sys.exit()


    def check_parameter_validity(self):
        """Check if various parameters are valid"""
        if self.number_of_phonons < self.output_trajectories_of_first:
            self.output_trajectories_of_first = self.number_of_phonons

        if self.number_of_timeframes <= self.number_of_stabilization_timeframes:
            logging.error("Parameter NUMBER_OF_STABILIZATION_TIMEFRAMES exceeds or equal to NUMBER_OF_TIMEFRAMES.\n" +
                          "Leave at least one timeframe after NUMBER_OF_STABILIZATION_TIMEFRAMES.\n" +
                          f"See the documentation at {WEBSITE}")
            sys.exit()

        for source in self.phonon_sources:
            if source.y > self.length:
                logging.error("Y coordinate of a source exceeded LENGHT")
                sys.exit()

            if source.y < 0:
                logging.error("Y coordinate of a source is negative.")
                sys.exit()

            if source.y - source.size_y / 2 < 0:
                logging.error("Source size along Y coordinate is too large")
                sys.exit()

            if abs(source.x) > self.width/2:
                logging.error("X coordinate of a source exceeds WIDTH")
                sys.exit()

            if abs(source.x + source.size_x / 2) > self.width/2:
                logging.error("Source size along X coordinate is too large")
                sys.exit()

            if abs(source.z) > self.thickness/2:
                logging.error("Z coordinate of a source exceeds THICKNESS")
                sys.exit()

            if abs(source.z + source.size_z / 2) > self.thickness/2:
                logging.error("Source size along Z coordinate is too large")
                sys.exit()

        if self.output_path_animation and self.number_of_timesteps > 5000:
            logging.warning("NUMBER_OF_TIMESTEPS is rather large for animation")

        if (self.cold_side_position_top and self.include_top_sidewall or
            self.hot_side_position_top and self.include_top_sidewall or
            self.cold_side_position_top and self.hot_side_position_top):
            logging.error("Top side is assigned multiple functions")
            sys.exit()

        if (self.cold_side_position_bottom and self.include_bottom_sidewall or
            self.hot_side_position_bottom and self.include_bottom_sidewall or
            self.cold_side_position_bottom and self.hot_side_position_bottom):
            logging.error("Bottom side is assigned multiple functions")
            sys.exit()

        if (self.cold_side_position_right and self.include_right_sidewall or
            self.hot_side_position_right and self.include_right_sidewall or
            self.cold_side_position_right and self.hot_side_position_right):
            logging.error("Right side is assigned multiple functions")
            sys.exit()

        if (self.cold_side_position_left and self.include_left_sidewall or
            self.hot_side_position_left and self.include_left_sidewall or
            self.cold_side_position_left and self.hot_side_position_left):
            logging.error("Left side is assigned multiple functions")
            sys.exit()


    def check_depricated_parameters(self):
        """Check for deprecated parameters and warn about them"""

        if 'COLD_SIDE_POSITION' in globals():
            logging.error("Parameter COLD_SIDE_POSITION is deprecated.\n" +
                          "Use specific boolean parameters like COLD_SIDE_POSITION_TOP = True\n" +
                          f"See the documentation at {WEBSITE}")
            sys.exit()

        if 'SPECIFIC_HEAT_CAPACITY' in globals():
            logging.warning("Parameter SPECIFIC_HEAT_CAPACITY is deprecated. Heat capacity is set from material.")

        if any([
            'HOT_SIDE_POSITION' in globals(),
            'HOT_SIDE_X' in globals(),
            'HOT_SIDE_Y' in globals(),
            'HOT_SIDE_WIDTH_X' in globals(),
            'HOT_SIDE_WIDTH_Y' in globals(),
            'HOT_SIDE_ANGLE_DISTRIBUTION' in globals(),
            'PHONON_SOURCE_X' in globals(),
            'PHONON_SOURCE_Y' in globals(),
            'PHONON_SOURCE_WIDTH_X' in globals(),
            'PHONON_SOURCE_WIDTH_Y' in globals(),
            'PHONON_SOURCE_ANGLE_DISTRIBUTION' in globals(),
            ]):
            logging.error("Parameters related to HOT_SIDE_... or PHONON_SOURCE_... were deprecated.\n" +
                          "Phonon source should be defined through the PHONON_SOURCES variable.\n" +
                          f"See the documentation at {WEBSITE}")
            sys.exit()

cf = Config()
cf.convert_to_enums()
cf.check_parameter_validity()
cf.check_depricated_parameters()
