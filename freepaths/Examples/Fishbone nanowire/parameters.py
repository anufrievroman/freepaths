"""Config file to simulate a fishbone nanowire in Si at 4K"""

import numpy as np
from options import *

# General parameters:
OUTPUT_FOLDER_NAME             = 'Fishbone nanowire'
NUMBER_OF_PHONONS              = 300
NUMBER_OF_TIMESTEPS            = 30000
NUMBER_OF_NODES                = 400
TIMESTEP                       = 1e-12
T                              = 4.0
PLOTS_IN_TERMINAL              = False
OUTPUT_SCATTERING_MAP          = False
OUTPUT_RAW_THERMAL_MAP         = True
OUTPUT_TRAJECTORIES_OF_FIRST   = 10
NUMBER_OF_LENGTH_SEGMENTS      = 10
HOT_SIDE_ANGLE_DISTRIBUTION    = Distributions.RANDOM


# Map & profiles parameters:
NUMBER_OF_PIXELS_X             = 100
NUMBER_OF_PIXELS_Y             = 100
NUMBER_OF_TIMEFRAMES           = 6


# Material parameters:
MEDIA                          = Materials.SILICON
SPECIFIC_HEAT_CAPACITY         = 0.0176  # [J/kg/K] for Si at 4 K (NOT CORRECT)
#SPECIFIC_HEAT_CAPACITY        = 714     # [J/kg/K] for Si at 300 K


# Internal scattering:
INCLUDE_INTERNAL_SCATTERING    = True
USE_GRAY_APPROXIMATION_MFP     = False
GRAY_APPROXIMATION_MFP         = 200e-9


# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 500e-9
LENGTH                         = 1500e-9


# Hot and cold sides [m]:
FREQUENCY_DETECTOR_SIZE        = WIDTH
COLD_SIZE_POSITION             = Positions.TOP
HOT_SIZE_X                     = 0
HOT_SIZE_WIDTH                 = 100e-9


# Roughness [m]:
SIDE_WALL_ROUGHNESS            = 2e-9
HOLE_ROUGHNESS                 = 2e-9
PILLAR_ROUGHNESS               = 2e-9
TOP_ROUGHNESS                  = 0.2e-9
BOTTOM_ROUGHNESS               = 0.2e-9
PILLAR_TOP_ROUGHNESS           = 0.2e-9


# Hole array parameters [m]:
INCLUDE_HOLES                  = True
CIRCULAR_HOLE_DIAMETER         = None
RECTANGULAR_HOLE_SIDE_X        = 300e-9
RECTANGULAR_HOLE_SIDE_Y        = 150e-9
PERIOD_X                       = WIDTH
PERIOD_Y                       = 300e-9

FIRST_HOLE_COORDINATE = 0
NUMBER_OF_PERIODS_X = 2
NUMBER_OF_PERIODS_Y = 6
HOLE_COORDINATES = np.zeros((NUMBER_OF_PERIODS_X * NUMBER_OF_PERIODS_Y, 3))
HOLE_SHAPES = [Shapes.RECTANGLE for x in range(HOLE_COORDINATES.shape[0])]
hole_number = 0
for i in range(NUMBER_OF_PERIODS_Y):
    for j in range(NUMBER_OF_PERIODS_X):
        HOLE_COORDINATES[hole_number, 0] = -(NUMBER_OF_PERIODS_X - 1) * PERIOD_X / 2 + j * PERIOD_X
        HOLE_COORDINATES[hole_number, 1] = FIRST_HOLE_COORDINATE + i * PERIOD_Y
        hole_number += 1


# Pillar array parameters [m]
INCLUDE_PILLARS = False
