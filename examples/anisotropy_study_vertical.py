"""Config file to simulate a membrane with a staggered lattice of rectangular slits
Here we impose thermal gradient in vertical direction"""

import numpy as np


# General parameters:
OUTPUT_FOLDER_NAME             = 'Anisotropy study vertical'
NUMBER_OF_PHONONS              = 1000
NUMBER_OF_TIMESTEPS            = 30000
NUMBER_OF_NODES                = 400
TIMESTEP                       = 1e-12
T                              = 4.0

# Map & profiles parameters:
NUMBER_OF_PIXELS_X             = 100
NUMBER_OF_PIXELS_Y             = 100
NUMBER_OF_TIMEFRAMES           = 6

# Material parameters:
MEDIA                          = 'Si'
SPECIFIC_HEAT_CAPACITY         = 0.0176  # [J/kg/K] for Si at 4 K


# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 1300e-9
LENGTH                         = 1300e-9

# STRUCTURE:

# |    C    |   
# |         |  ^
# |         |  |
# |         |  |
# |    H    | 

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

# Phonon source:
PHONON_SOURCES = [Source(x=0, y=0, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="random_up")]

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
