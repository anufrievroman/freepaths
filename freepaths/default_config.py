"""Default config file"""

import numpy as np

# General parameters:
OUTPUT_FOLDER_NAME             = "Si nanowire at 300 K"
NUMBER_OF_PHONONS              = 500
NUMBER_OF_TIMESTEPS            = 60000
NUMBER_OF_NODES                = 400
TIMESTEP                       = 2e-12
T                              = 300
PLOTS_IN_TERMINAL              = False
OUTPUT_SCATTERING_MAP          = False
OUTPUT_RAW_THERMAL_MAP         = True
OUTPUT_TRAJECTORIES_OF_FIRST   = 10
NUMBER_OF_LENGTH_SEGMENTS      = 10
HOT_SIDE_ANGLE_DISTRIBUTION    = "random"

# Animation:
OUTPUT_PATH_ANIMATION          = False
PHONONS_IN_ANIMATION           = 3

# Map & profiles parameters:
NUMBER_OF_PIXELS_X             = 100
NUMBER_OF_PIXELS_Y             = 100
NUMBER_OF_TIMEFRAMES           = 5

# Material parameters:
MEDIA                          = "Si"
SPECIFIC_HEAT_CAPACITY         = 714  # [J/kg/K] for Si at 300 K

# Internal scattering:
INCLUDE_INTERNAL_SCATTERING    = True
USE_GRAY_APPROXIMATION_MFP     = False
GRAY_APPROXIMATION_MFP         = None

# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 500e-9
LENGTH                         = 2200e-9

# Hot and cold sides [m]:
FREQUENCY_DETECTOR_SIZE        = WIDTH
COLD_SIDE_POSITION             = 'top'
HOT_SIDE_X                     = 0
HOT_SIDE_WIDTH_X               = WIDTH
HOT_SIDE_Y                     = 0
HOT_SIDE_WIDTH_Y               = 0

# Roughness [m]:
SIDE_WALL_ROUGHNESS            = 2e-9
HOLE_ROUGHNESS                 = 2e-9
PILLAR_ROUGHNESS               = 2e-9
TOP_ROUGHNESS                  = 0.2e-9
BOTTOM_ROUGHNESS               = 0.2e-9
PILLAR_TOP_ROUGHNESS           = 0.2e-9

# Parabolic boundaries:
INCLUDE_TOP_PARABOLA           = False
TOP_PARABOLA_TIP               = 1000e-9
TOP_PARABOLA_FOCUS             = 100e-9
INCLUDE_BOTTOM_PARABOLA        = False
BOTTOM_PARABOLA_TIP            = 1000e-9
BOTTOM_PARABOLA_FOCUS          = 100e-9

# Hole array parameters [m]:
INCLUDE_HOLES                  = False
CIRCULAR_HOLE_DIAMETER         = 200e-9
RECTANGULAR_HOLE_SIDE_X        = 200e-9
RECTANGULAR_HOLE_SIDE_Y        = 200e-9
PERIOD_X                       = 300e-9
PERIOD_Y                       = 300e-9

# Lattice of holes:
HOLE_COORDINATES = np.zeros((1, 3))
HOLE_SHAPES = ["circle" for x in range(1)]

# Pillar array parameters [m]
INCLUDE_PILLARS                = False
PILLAR_HEIGHT                  = 100e-9
PILLAR_WALL_ANGLE              = 0.0
PILLAR_COORDINATES = np.zeros((1, 3))
