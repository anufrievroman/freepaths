# GENERAL PARAMETERS
output_folder_name='diode_4K_up'
number_of_phonons=5000    
number_of_phonons_in_a_group=50	    # To reduce the memory usage, phonons are simulated in small groups
number_of_timesteps=200000	    # Maximum number of timesteps
number_of_nodes=400		    # Resolution of distribution plots
timestep=1.0e-12                    # [s] Duration of one timestep
T=4.0                              # [K] Temperature of the system
output_in_terminal=False	    # Prevents showing plots directly in the terminal
output_scattering_map=False	    # Prevents calculating and outputing heavy scattering map
output_raw_thermal_map=False	    # Prevents outputing heavy raw data for thermal map

# SYSTEM DIMENSIONS [m]
width=4*800e-9
length=4*800e-9 
thickness=145e-9

# ROUGHNESS [m]
side_wall_roughness=0.01e-9
hole_roughness=2.0e-9
pillar_roughness=2.0e-9
top_roughness=0.02e-9
bottom_roughness=0.02e-9
pillar_top_roughness=2.0e-9

# HOLES AND PILLARS PARAMETERS [m]
holes='yes'                                                             
hole_lattice_type='staggered_triangles_up'
pillars='no'
pillar_lattice_type='square'
circular_hole_diameter=20e-9
rectangular_hole_side_x=440e-9
rectangular_hole_side_y=720e-9
pillar_height=30e-9
from math import pi
pillar_wall_angle=pi/2.0
period_x=800e-9
period_y=800e-9

# ENERGY MAP PARAMETERS
number_of_pixels_x=50	
number_of_pixels_y=50
number_of_timeframes=10

# MATERIAL PARAMETERS
specific_heat_capacity=714  #[J/kg/K] It'  0.0176 at 4K 
material_density=2330	    # [kg/m^3]
