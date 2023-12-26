"""Config file to simulate a fishbone nanowire in Si at 4K"""

import numpy as np

# General parameters:
OUTPUT_FOLDER_NAME             = 'Fishbone nanowire'
NUMBER_OF_PHONONS              = 300
NUMBER_OF_TIMESTEPS            = 30000
NUMBER_OF_NODES                = 400
TIMESTEP                       = 1e-12
T                              = 4.0
OUTPUT_TRAJECTORIES_OF_FIRST   = 30
NUMBER_OF_LENGTH_SEGMENTS      = 10


# Map & profiles parameters:
NUMBER_OF_PIXELS_X             = 100
NUMBER_OF_PIXELS_Y             = 100
NUMBER_OF_TIMEFRAMES           = 6


# Material parameters:
MEDIA                          = "Si"
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


# Phonon source:
PHONON_SOURCES                 = [Source(x=0, y=0, z=0, size_x=100e-9,  size_y=0, size_z=THICKNESS, angle_distribution="random_up")]


# Roughness [m]:
SIDE_WALL_ROUGHNESS            = 2e-9
HOLE_ROUGHNESS                 = 2e-9
TOP_ROUGHNESS                  = 0.2e-9
BOTTOM_ROUGHNESS               = 0.2e-9


# Hole array parameters [m]:
HOLES = []

PERIOD_X                       = WIDTH
PERIOD_Y                       = 300e-9
NUMBER_OF_PERIODS_X = 2
NUMBER_OF_PERIODS_Y = 6
for i in range(NUMBER_OF_PERIODS_Y):
    for j in range(NUMBER_OF_PERIODS_X):
        x_coor = -(NUMBER_OF_PERIODS_X - 1) * PERIOD_X / 2 + j * PERIOD_X
        y_coor = i * PERIOD_Y
        HOLES.append(RectangularHole(x=x_coor, y=y_coor, size_x=300e-9, size_y=150e-9))
