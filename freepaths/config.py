"""Module that reads the user input file, provides default values, and converts the variables into enums"""

import sys
import argparse

from freepaths.options import Materials, Distributions, Positions

# Import a default input file:
from freepaths.default_config import *


# Parse user arguments:
WEBSITE = 'https://anufrievroman.gitbook.io/freepaths'
parser = argparse.ArgumentParser(prog='FreePATHS', description='Monte Carlo simulator',
                                 epilog=f'For more information, visit: {WEBSITE}')
parser.add_argument('input_file', nargs='?', default=None, help='The input file')
parser.add_argument("-s", "--sampling", help="Run in MFP sampling mode", action="store_true")
args = parser.parse_args()


# If a file is provided, overwrite the default values:
if args.input_file:
    exec(open(args.input_file, encoding='utf-8').read(), globals())
else:
    print("You didn't provide any input file, so let's run a demo simulation!\n")


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
        self.output_structure_color = OUTPUT_STRUCTURE_COLOR
        self.number_of_length_segments = NUMBER_OF_LENGTH_SEGMENTS
        self.hot_side_angle_distribution = HOT_SIDE_ANGLE_DISTRIBUTION

        # Animation:
        self.output_path_animation = OUTPUT_PATH_ANIMATION
        self.output_animation_fps = OUTPUT_ANIMATION_FPS

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
        self.include_right_sidewall = INCLUDE_RIGHT_SIDEWALL
        self.include_left_sidewall = INCLUDE_LEFT_SIDEWALL
        self.include_top_sidewall = INCLUDE_TOP_SIDEWALL
        self.include_bottom_sidewall = INCLUDE_BOTTOM_SIDEWALL

        # Hot and cold sides:
        self.frequency_detector_size = FREQUENCY_DETECTOR_SIZE
        self.cold_side_position = COLD_SIDE_POSITION
        self.hot_side_position = HOT_SIDE_POSITION
        self.hot_side_x = HOT_SIDE_X
        self.hot_side_y = HOT_SIDE_Y
        self.hot_side_width_x = HOT_SIDE_WIDTH_X
        self.hot_side_width_y = HOT_SIDE_WIDTH_Y

        # Roughness:
        self.side_wall_roughness = SIDE_WALL_ROUGHNESS
        self.hole_roughness = HOLE_ROUGHNESS
        self.pillar_roughness = PILLAR_ROUGHNESS
        self.top_roughness = TOP_ROUGHNESS
        self.bottom_roughness = BOTTOM_ROUGHNESS
        self.pillar_top_roughness = PILLAR_TOP_ROUGHNESS

        # Parabolic boundary:
        self.include_top_parabola = INCLUDE_TOP_PARABOLA
        self.top_parabola_tip = TOP_PARABOLA_TIP
        self.top_parabola_focus = TOP_PARABOLA_FOCUS
        self.include_bottom_parabola = INCLUDE_BOTTOM_PARABOLA
        self.bottom_parabola_tip = BOTTOM_PARABOLA_TIP
        self.bottom_parabola_focus = BOTTOM_PARABOLA_FOCUS

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

        # Distributions:
        valid_distributions =[member.name.lower() for member in Distributions]
        if self.hot_side_angle_distribution in valid_distributions:
            self.hot_side_angle_distribution = Distributions[self.hot_side_angle_distribution.upper()]
        else:
            print("ERROR: Parameter HOT_SIDE_ANGLE_DISTRIBUTION is not set correctly.")
            print("HOT_SIDE_ANGLE_DISTRIBUTION should be one of the following:")
            print(*valid_distributions, sep = ", ")
            sys.exit()

        # Materials:
        valid_materials = [member.name for member in Materials]
        if self.media in valid_materials:
            self.media = Materials[self.media]
        else:
            print(f"ERROR: Material {self.media} is not in the database.")
            print("MEDIA should be one of the following:")
            print(*valid_materials, sep = ", ")
            sys.exit()

        # Positions:
        valid_positions = [member.name.lower() for member in Positions]
        if self.cold_side_position in valid_positions:
            self.cold_side_position = Positions[self.cold_side_position.upper()]
        else:
            print("ERROR: Parameter COLD_SIDE_POSITION is not set correctly.")
            print("COLD_SIDE_POSITION should be one of the following:")
            print(*valid_positions, sep = ", ")
            sys.exit()

        if self.hot_side_position in valid_positions:
            self.hot_side_position = Positions[self.hot_side_position.upper()]
        else:
            print("ERROR: Parameter HOT_SIDE_POSITION is not set correctly.")
            print("HOT_SIDE_POSITION should be one of the following:")
            print(*valid_positions, sep = ", ")
            sys.exit()


    def check_parameter_validity(self):
        """Check if various parameters are valid"""
        if self.number_of_phonons < self.output_trajectories_of_first:
            self.output_trajectories_of_first = self.number_of_phonons
            print("WARNING: Parameter OUTPUT_TRAJECTORIES_OF_FIRST exceeded NUMBER_OF_PHONONS.\n")

        if self.hot_side_y > self.length:
            self.hot_side_y = self.length
            print("WARNING: Parameter HOT_SIDE_Y exceeded LENGHT.\n")

        if self.hot_side_y < 0:
            self.hot_side_y = 0
            print("WARNING: Parameter HOT_SIDE_Y was negative.\n")

        if self.hot_side_y - self.hot_side_width_y / 2 < 0:
            self.hot_side_width_y = self.hot_side_y * 2
            print("WARNING: Parameter HOT_SIDE_WIDTH_Y was too large.\n")

        if self.hot_side_x > self.width/2:
            self.hot_side_x = 0
            print("WARNING: Parameter HOT_SIDE_X was larger than WIDTH.\n")

        if self.hot_side_width_x > self.width:
            self.hot_side_width_x = self.width
            print("WARNING: Parameter HOT_SIDE_WIDTH_X exceeds WIDTH.\n")

        if self.cold_side_position == self.hot_side_position:
            print(f"ERROR: Hot and cold sides are set at {self.cold_side_position}.")
            sys.exit()

        if self.output_path_animation and self.number_of_timesteps > 5000:
            print("WARNING: NUMBER_OF_TIMESTEPS is rather large for animation.\n")


cf = Config()
cf.convert_to_enums()
cf.check_parameter_validity()
