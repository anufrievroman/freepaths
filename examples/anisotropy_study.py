"""
Config file to simulate a membrane with a staggered lattice of rectangular slits
Here we impose thermal gradient in vertical or horizonta direction
"""

import numpy as np


# General parameters:
OUTPUT_FOLDER_NAME             = 'Anisotropy study'
NUMBER_OF_PARTICLES            = 2000
MEDIA                          = 'Si'
T                              = 300.0

# Multiprocessing
NUMBER_OF_PROCESSES = 10

# Simulation time parameters:
TIMESTEP                       = 1.0e-12
NUMBER_OF_TIMESTEPS            = 100000
NUMBER_OF_VIRTUAL_TIMESTEPS    = NUMBER_OF_TIMESTEPS*4
NUMBER_OF_STABILIZATION_TIMEFRAMES = 5
NUMBER_OF_TIMEFRAMES = 8

# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 1300e-9
LENGTH                         = 1300e-9

# Map & profiles parameters:
pixel_size = 20e-9
NUMBER_OF_PIXELS_X             = int(WIDTH / pixel_size)
NUMBER_OF_PIXELS_Y             = int(LENGTH / pixel_size)
IGNORE_FAULTY_PARTICLES        = False


# VERTICAL GRADIENT:

# Walls:
INCLUDE_RIGHT_SIDEWALL         = True
INCLUDE_LEFT_SIDEWALL          = True
INCLUDE_TOP_SIDEWALL           = False
INCLUDE_BOTTOM_SIDEWALL        = False

# Hot and cold sides [m]:
COLD_SIDE_POSITION_RIGHT       = False
COLD_SIDE_POSITION_TOP         = True
HOT_SIDE_POSITION_LEFT         = False
HOT_SIDE_POSITION_BOTTOM       = True

# Particle source:
PARTICLE_SOURCES = [Source(x=0, y=0, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="random", angle=0)]


# HORIZONTAL GRADIENT:

# INCLUDE_RIGHT_SIDEWALL         = False
# INCLUDE_LEFT_SIDEWALL          = False
# INCLUDE_TOP_SIDEWALL           = True
# INCLUDE_BOTTOM_SIDEWALL        = True

# Hot and cold sides [m]:
# COLD_SIDE_POSITION_TOP              = False
# COLD_SIDE_POSITION_RIGHT            = True
# HOT_SIDE_POSITION_LEFT              = True
# HOT_SIDE_POSITION_BOTTOM            = False

# Particle source:
PARTICLE_SOURCES = [Source(x=-WIDTH/2, y=LENGTH/2, z=0, size_x=0,  size_y=LENGTH, size_z=THICKNESS, angle_distribution="random", angle=np.pi/2)]



HOLES = []

# Staggered attice of holes:
period_x = 300e-9
period_y = 300e-9
number_of_periods_x = 4
number_of_periods_y = 4
for i in range(number_of_periods_y):
    for j in range(number_of_periods_x):
        x = -(number_of_periods_x - 1) * period_x / 2 + j * period_x
        y = 200e-9 + i * period_y
        HOLES.append(RectangularHole(x, y, size_x=200e-9, size_y=100e-9))
