from math import sin, pi, cos

# General parameters:
OUTPUT_FOLDER_NAME             = 'Test'
NUMBER_OF_PHONONS              = 100
NUMBER_OF_TIMESTEPS            = 30000
NUMBER_OF_NODES                = 400
TIMESTEP                       = 1.0e-12
T                              = 4.0
PLOTS_IN_TERMINAL              = False
OUTPUT_SCATTERING_MAP          = True
OUTPUT_RAW_THERMAL_MAP         = True
OUTPUT_TRAJECTORIES_OF_FIRST   = 100
FREQUENCY_DETECTOR_SIZE        = 3000e-9
COLD_SIZE_POSITION             = 'top'
NUMBER_OF_LENGTH_SEGMENTS      = 10
HOT_SIDE_ANGLE_DISTRIBUTION    = "random"

# Internal scattering:
INCLUDE_INTERNAL_SCATTERING    = True
USE_GRAY_APPROXIMATION_MFP     = False
GRAY_APPROXIMATION_MFP         = 200e-9

# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 300e-9
LENGTH                         = 1000e-9

# Roughness [m]:
SIDE_WALL_ROUGHNESS            = 2e-9
HOLE_ROUGHNESS                 = 2e-9
PILLAR_ROUGHNESS               = 2e-9
TOP_ROUGHNESS                  = 0.2e-9
BOTTOM_ROUGHNESS               = 0.2e-9
PILLAR_TOP_ROUGHNESS           = 0.2e-9

# Holes and pillars parameters [m]:
INCLUDE_HOLES                  = True
HOLE_LATTICE_TYPE              = 'square'
INCLUDE_PILLARS                = False
PILLAR_LATTICE_TYPE            = 'square'
CIRCULAR_HOLE_DIAMETER         = 150e-9
RECTANGULAR_HOLE_SIDE_X        = 2 * THICKNESS * sin(pi / 6)
RECTANGULAR_HOLE_SIDE_Y        = THICKNESS * cos(pi / 6)
PILLAR_HEIGHT                  = 30e-9
PILLAR_WALL_ANGLE              = pi / 2.0
PERIOD_X                       = 300e-9
PERIOD_Y                       = 300e-9

# Map & profiles parameters:
NUMBER_OF_PIXELS_X             = 100
NUMBER_OF_PIXELS_Y             = 100
NUMBER_OF_TIMEFRAMES           = 6

# Material parameters:
MEDIA                          = 'Si'
SPECIFIC_HEAT_CAPACITY         = 714  # [J/kg/K] for Si at 300 K
# specific_heat_capacity       = 606  # [J/kg/K] for SiC at 300 K
# specific_heat_capacity       = 0.0176  # [J/kg/K] for Si at 4 K
