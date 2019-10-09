from math import pi

# SIMULATION PARAMETERS
# You must create the output_folder before you start the simulation
output_folder_name='membrane'                                                       
number_of_phonons=2000    
number_of_phonons_in_a_group=100                                                 # To reduce the memory load, phonons are simulated in small groups
number_of_timesteps=20000
number_of_nodes=400                                                             # Resolution of distribution plots
timestep=0.5e-12                                                                # [s] Duration of one timestep
T=300.0                                                                         # [K] Temperature of the system
output_in_terminal=False

# SYSTEM DIMENSIONS [m]
width=300e-9
length=300e-9 #2*600e-9 - 3*155e-9
thickness=50e-9

# ROUGHNESS [m]
side_wall_roughness=2.0e-9
hole_roughness=1.0e-9
pillar_roughness=2.0e-9
top_roughness=0.02e-9
bottom_roughness=0.02e-9
pillar_top_roughness=2.0e-9

# SCATTER PARAMETERS [m]
holes='no'                                                             
hole_lattice_type='square'#'square'
pillars='no'
pillar_lattice_type='black_silicon'
circular_hole_diameter=20e-9#185e-9
rectangular_hole_side_x=500e-9#600e-9 - 155e-9
rectangular_hole_side_y=100e-9#600e-9 - 2*155e-9                                                        
pillar_height=30e-9
pillar_wall_angle=pi/2.0 #pi/3.0                                                          
period_x=50e-9
period_y=50e-9

# ENERGY MAP PARAMETERS
number_of_pixels_x=int(0.2*width*1e9)
number_of_pixels_y=int(0.2*length*1e9)
number_of_timeframes=20

# MATERIAL PARAMETERS
specific_heat_capacity=714 #0.0176 #714 #[J/kg/K]
material_density=2330 #[kg/m^3]
