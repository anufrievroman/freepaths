# GENERAL PARAMETERS
output_folder_name = 'Source n_200nm'      # This folder will be created 
number_of_phonons = 40000            # Total number of phonons that will be simulated
number_of_phonons_in_a_group = 100  # To reduce the memory usage, phonons are simulated in small groups
number_of_timesteps = 150000	    # Maximum number of timesteps
number_of_nodes = 400		    # Resolution of distribution plots
timestep = 1.0e-12                  # [s] Duration of one timestep
T = 4.0		                    # [K] Temperature of the system
output_in_terminal = False	    # Prevents showing plots directly in the terminal
output_scattering_map = False	    # Prevents calculating and outputing heavy scattering map
output_raw_thermal_map = False	    # Prevents outputing heavy raw data for thermal map
simulation_mode = 1

# INTERNAL SCATTERING
internal_scattering_on = True      
use_gray_approximation_mfp = False   # Do you want to use frequency indepenednt MFP specified below?
gray_approximation_mfp = 1.5e-6      # [m] Frequency independent internal MFP

# SYSTEM DIMENSIONS [m]
width = 4*300e-9
length = 12*300e-9 
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
circular_hole_diameter = 100e-9
rectangular_hole_side_x = 250e-9
rectangular_hole_side_y = 20e-9
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
