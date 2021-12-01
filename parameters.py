# GENERAL PARAMETERS
output_folder_name           = 'Test'
number_of_phonons            = 1000
number_of_phonons_in_a_group = 200
number_of_timesteps          = 40000
number_of_nodes              = 400
timestep                     = 1.0e-12
T                            = 4.0
output_in_terminal           = False
output_scattering_map        = False
output_raw_thermal_map       = True
frequency_detector_size      = 3000e-9
cold_side                    = 'top'
simulation_mode              = 1
number_of_length_segments    = 10
hot_side_angle_distribution  = "random"

# INTERNAL SCATTERING
internal_scattering_on       = True
use_gray_approximation_mfp   = False
gray_approximation_mfp       = 500e-9

# SYSTEM DIMENSIONS [m]
width                        = 1000e-9
length                       = 3000e-9
thickness                    = 100e-9

# ROUGHNESS [m]
side_wall_roughness          = 2.0e-9
hole_roughness               = 2.0e-9
pillar_roughness             = 2.0e-9
top_roughness                = 2.0e-9
bottom_roughness             = 2.0e-9
pillar_top_roughness         = 2.0e-9

# HOLES AND PILLARS PARAMETERS [m]
holes                        = 'yes'
hole_lattice_type            = 'square'
pillars                      = 'no'
pillar_lattice_type          = 'square'
circular_hole_diameter       = 500e-9
rectangular_hole_side_x      = 500e-9
rectangular_hole_side_y      = 500e-9
pillar_height                = 30e-9
import math
pillar_wall_angle            = math.pi/2.0
period_x                     = 300e-9
period_y                     = 1000e-9

# MAP & PROFILES PARAMETERS
number_of_pixels_x           = 150
number_of_pixels_y           = 150
number_of_timeframes         = 6

# MATERIAL PARAMETERS
material                     = 'Si'
specific_heat_capacity       = 714  # [J/kg/K] for Si at 300 K
# specific_heat_capacity       = 606  # [J/kg/K] for SiC at 300 K
# specific_heat_capacity       = 0.0176  # [J/kg/K] for Si at 4 K
