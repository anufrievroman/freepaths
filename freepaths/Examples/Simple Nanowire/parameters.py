"""Config file to simulate a simple nanowire"""

import numpy as np


# General parameters:
OUTPUT_FOLDER_NAME             = 'Nanowire'
NUMBER_OF_PHONONS              = 1000
NUMBER_OF_TIMESTEPS            = 30000
NUMBER_OF_NODES                = 400
TIMESTEP                       = 1.0e-12
T                              = 4.0
PLOTS_IN_TERMINAL              = False
OUTPUT_SCATTERING_MAP          = False
OUTPUT_RAW_THERMAL_MAP         = True
OUTPUT_TRAJECTORIES_OF_FIRST   = 50
NUMBER_OF_LENGTH_SEGMENTS      = 10
HOT_SIDE_ANGLE_DISTRIBUTION    = "random"


# Map & profiles parameters:
NUMBER_OF_PIXELS_X             = 20
NUMBER_OF_PIXELS_Y             = 100
NUMBER_OF_TIMEFRAMES           = 6


# Material parameters:
MEDIA                          = 'Si'
SPECIFIC_HEAT_CAPACITY         = 0.0176  # [J/kg/K] for Si at 4 K
#SPECIFIC_HEAT_CAPACITY        = 714     # [J/kg/K] for Si at 300 K


# Internal scattering:
INCLUDE_INTERNAL_SCATTERING    = True
USE_GRAY_APPROXIMATION_MFP     = False
GRAY_APPROXIMATION_MFP         = None


# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 200e-9
LENGTH                         = 2000e-9

# Hot and cold sides [m]:
FREQUENCY_DETECTOR_SIZE        = WIDTH
COLD_SIZE_POSITION             = 'top'
HOT_SIZE_X                     = 0
HOT_SIZE_WIDTH                 = WIDTH

# Roughness [m]:
SIDE_WALL_ROUGHNESS            = 2e-9
HOLE_ROUGHNESS                 = 2e-9
PILLAR_ROUGHNESS               = 2e-9
TOP_ROUGHNESS                  = 0.2e-9
BOTTOM_ROUGHNESS               = 0.2e-9
PILLAR_TOP_ROUGHNESS           = 0.2e-9


# Hole array parameters [m]:
INCLUDE_HOLES                  = False

# Pillar array parameters [m]
INCLUDE_PILLARS                = False
