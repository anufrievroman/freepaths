# General parameters:
output_folder_name           = 'Example'
number_of_phonons            = 1000
number_of_phonons_in_a_group = 100
number_of_timesteps          = 30000
number_of_nodes              = 400
timestep                     = 1.0e-12
T                            = 4.0
output_in_terminal           = False
output_scattering_map        = False
output_raw_thermal_map       = False
frequency_detector_size      = 5000e-9
cold_side                    = 'top'
simulation_mode              = 1
number_of_length_segments    = 10
hot_side_angle_distribution  = 'random'

# Internal scattering:
internal_scattering_on       = True
use_gray_approximation_mfp   = False
gray_approximation_mfp       = 1600e-9

# System dimensions [m]:
from math import sin, pi, cos
thickness                    = 1000e-9
width                        = 600e-9
length                       = 1500e-9

# Roughness [m]:
side_wall_roughness          = 2e-9
hole_roughness               = 2e-9
pillar_roughness             = 2e-9
top_roughness                = 0.2e-9
bottom_roughness             = 0.2e-9
pillar_top_roughness         = 0.2e-9

# Holes and pillars parameters [m]:
holes                        = 'yes'
hole_lattice_type            = 'square'
pillars                      = 'no'
pillar_lattice_type          = 'square'
circular_hole_diameter       = 200e-9
rectangular_hole_side_x      = 200e-9
rectangular_hole_side_y      = 200e-9
pillar_height                = 30e-9
import math
pillar_wall_angle            = math.pi/2.0
period_x                     = 300e-9
period_y                     = 300e-9

# Map & profiles parameters:
number_of_pixels_x           = 100
number_of_pixels_y           = 100
number_of_timeframes         = 6

# Material parameters:
material                     = 'Si'
# specific_heat_capacity       = 714  # [J/kg/K] for Si at 300 K
# specific_heat_capacity       = 606  # [J/kg/K] for SiC at 300 K
specific_heat_capacity       = 0.0176  # [J/kg/K] for Si at 4 K
