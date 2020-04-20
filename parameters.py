# GENERAL PARAMETERS
output_folder_name = 'time test'
number_of_phonons = 20
number_of_phonons_in_a_group = 20
number_of_timesteps = 150000
number_of_nodes = 400
timestep = 1.0e-12      
T = 4.0		         
output_in_terminal = False
output_scattering_map = False
output_raw_thermal_map = True
calculate_mfp_spectrum = False
calculate_travel_times = False
cold_side_on_top = True
simulation_mode = 1

# INTERNAL SCATTERING
internal_scattering_on = True      
use_gray_approximation_mfp = False
gray_approximation_mfp = 1.5e-6

# SYSTEM DIMENSIONS [m]
width = 1200e-9
length = 3200e-9 
thickness = 145e-9

# ROUGHNESS [m]
side_wall_roughness = 2.0e-12
hole_roughness = 2.0e-9
pillar_roughness = 2.0e-9
top_roughness = 0.2e-9
bottom_roughness = 0.2e-9
pillar_top_roughness = 2.0e-9

# HOLES AND PILLARS PARAMETERS [m]
holes = 'yes'                                                             
hole_lattice_type = 'source'
pillars = 'no'
pillar_lattice_type = 'square'
circular_hole_diameter = 50e-9
rectangular_hole_side_x = 2400e-9
rectangular_hole_side_y = 1200e-9
pillar_height = 30e-9
from math import pi
pillar_wall_angle = pi/2.0
period_x = 300e-9
period_y = 300e-9

# ENERGY MAP PARAMETERS
number_of_pixels_x = 300
number_of_pixels_y = 300
number_of_timeframes = 10

# MATERIAL PARAMETERS
specific_heat_capacity = 714            #[J/kg/K] It's  0.0176 at 4K 
material_density = 2330	                # [kg/m^3]
