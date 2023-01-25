"""Config file to simulate a phononic crystal with square lattice of holes"""

import numpy as np


# General parameters:
OUTPUT_FOLDER_NAME             = 'Phononic crystal square'
NUMBER_OF_PHONONS              = 300
NUMBER_OF_TIMESTEPS            = 30000
NUMBER_OF_NODES                = 400
TIMESTEP                       = 1e-12
T                              = 4.0
PLOTS_IN_TERMINAL              = False
OUTPUT_SCATTERING_MAP          = False
OUTPUT_RAW_THERMAL_MAP         = True
OUTPUT_TRAJECTORIES_OF_FIRST   = 30
NUMBER_OF_LENGTH_SEGMENTS      = 10
HOT_SIDE_ANGLE_DISTRIBUTION    = 'random'


# Map & profiles parameters:
NUMBER_OF_PIXELS_X             = 100
NUMBER_OF_PIXELS_Y             = 100
NUMBER_OF_TIMEFRAMES           = 6


# Material parameters:
MEDIA                          = 'Si'
SPECIFIC_HEAT_CAPACITY         = 0.0176  # [J/kg/K] for Si at 4 K
#SPECIFIC_HEAT_CAPACITY        = 714  # [J/kg/K] for Si at 300 K
# SPECIFIC_HEAT_CAPACITY       = 606  # [J/kg/K] for SiC at 300 K


# Internal scattering:
INCLUDE_INTERNAL_SCATTERING    = True
USE_GRAY_APPROXIMATION_MFP     = False
GRAY_APPROXIMATION_MFP         = None


# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 1200e-9
LENGTH                         = 2200e-9

# Hot and cold sides [m]:
FREQUENCY_DETECTOR_SIZE        = WIDTH
COLD_SIZE_POSITION             = 'top'
HOT_SIZE_X                     = 0
HOT_SIZE_WIDTH                 = WIDTH


# Roughness [m]:
SIDE_WALL_ROUGHNESS            = 2e-9
HOLE_ROUGHNESS                 = 2e-9
TOP_ROUGHNESS                  = 0.2e-9
BOTTOM_ROUGHNESS               = 0.2e-9


# Hole array parameters [m]:
INCLUDE_HOLES                  = True
CIRCULAR_HOLE_DIAMETER         = 200e-9
RECTANGULAR_HOLE_SIDE_X        = None
RECTANGULAR_HOLE_SIDE_Y        = None
PERIOD_X                       = 300e-9
PERIOD_Y                       = 300e-9


# Lattice of holes:
FIRST_HOLE_COORDINATE = 300e-9
NUMBER_OF_PERIODS_X = 5
NUMBER_OF_PERIODS_Y = 6
HOLE_COORDINATES = np.zeros((NUMBER_OF_PERIODS_X * NUMBER_OF_PERIODS_Y, 3))
HOLE_SHAPES = ['circle' for x in range(HOLE_COORDINATES.shape[0])]
hole_number = 0
for i in range(NUMBER_OF_PERIODS_Y):
    for j in range(NUMBER_OF_PERIODS_X):
        HOLE_COORDINATES[hole_number, 0] = -(NUMBER_OF_PERIODS_X - 1) * PERIOD_X / 2 + j * PERIOD_X
        HOLE_COORDINATES[hole_number, 1] = FIRST_HOLE_COORDINATE + i * PERIOD_Y
        HOLE_COORDINATES[hole_number, 2] = 0
        hole_number += 1
