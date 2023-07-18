"""Config file to simulate a membrane with a staggered lattice of rectangular slits"""

import numpy as np


# General parameters:
OUTPUT_FOLDER_NAME             = 'Anisotropy study horizontal'
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

#      ---->
#-------------------
#
# H               C
#
#-------------------

INCLUDE_RIGHT_SIDEWALL         = False
INCLUDE_LEFT_SIDEWALL          = False
INCLUDE_TOP_SIDEWALL           = True
INCLUDE_BOTTOM_SIDEWALL        = True

# Hot and cold sides [m]:
COLD_SIDE_POSITION_TOP              = False
COLD_SIDE_POSITION_RIGHT            = True
HOT_SIDE_POSITION_LEFT              = True
HOT_SIDE_POSITION_BOTTOM            = False

# Phonon source:
PHONON_SOURCES = [Source(x=-WIDTH/2, y=LENGTH/2, z=0, size_x=0,  size_y=LENGTH, size_z=THICKNESS, angle_distribution="random_right")]

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
