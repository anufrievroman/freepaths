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
COLD_SIDE_POSITION             = 'right'
HOT_SIDE_POSITION              = 'left'
HOT_SIDE_X                     = -WIDTH/2
HOT_SIDE_Y                     = LENGTH/2
HOT_SIDE_WIDTH_X               = 0
HOT_SIDE_WIDTH_Y               = LENGTH
HOT_SIDE_ANGLE_DISTRIBUTION    = 'random_right'


# Hole array parameters [m]:
INCLUDE_HOLES                  = True
CIRCULAR_HOLE_DIAMETER         = 0
RECTANGULAR_HOLE_SIDE_X        = 200e-9
RECTANGULAR_HOLE_SIDE_Y        = 100e-9
PERIOD_X                       = 300e-9
PERIOD_Y                       = 300e-9


# Staggered attice of holes:
FIRST_HOLE_COORDINATE = 200e-9
NUMBER_OF_PERIODS_X = 4
NUMBER_OF_PERIODS_Y = 4
HOLE_COORDINATES = np.zeros((NUMBER_OF_PERIODS_X * NUMBER_OF_PERIODS_Y, 3))
HOLE_SHAPES = ['rectangle' for x in range(HOLE_COORDINATES.shape[0])]
hole_number = 0
for i in range(NUMBER_OF_PERIODS_Y):
    for j in range(NUMBER_OF_PERIODS_X):
        HOLE_COORDINATES[hole_number, 0] = -(NUMBER_OF_PERIODS_X - 1) * PERIOD_X / 2 + j * PERIOD_X
        HOLE_COORDINATES[hole_number, 1] = FIRST_HOLE_COORDINATE + i * PERIOD_Y
        HOLE_COORDINATES[hole_number, 2] = 0
        hole_number += 1
