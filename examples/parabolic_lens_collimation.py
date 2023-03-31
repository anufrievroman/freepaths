"""Config file to simulate a structure with parabolic surface at the bottom
which acts like a collimator for phonon flux emitted from the hot spot of 10x10nm"""

# General parameters:
OUTPUT_FOLDER_NAME             = 'Parabolic lens collimation'
NUMBER_OF_PHONONS              = 2000
NUMBER_OF_TIMESTEPS            = 30000
NUMBER_OF_NODES                = 400
TIMESTEP                       = 1.0e-12
T                              = 4.0
PLOTS_IN_TERMINAL              = False
OUTPUT_SCATTERING_MAP          = False
OUTPUT_RAW_THERMAL_MAP         = True
OUTPUT_TRAJECTORIES_OF_FIRST   = 30
NUMBER_OF_LENGTH_SEGMENTS      = 10


# Map & profiles parameters:
NUMBER_OF_PIXELS_X             = 100
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
WIDTH                          = 2000e-9
LENGTH                         = 2000e-9

# Hot and cold sides [m]:
FREQUENCY_DETECTOR_SIZE        = WIDTH
COLD_SIDE_POSITION             = 'top'
HOT_SIDE_X                     = 0
HOT_SIDE_Y                     = 300e-9
HOT_SIDE_WIDTH_X               = 10e-9
HOT_SIDE_WIDTH_Y               = 10e-9
HOT_SIDE_ANGLE_DISTRIBUTION    = 'uniform'

# Roughness [m]:
SIDE_WALL_ROUGHNESS            = 2e-9
TOP_ROUGHNESS                  = 0.2e-9
BOTTOM_ROUGHNESS               = 0.2e-9

# Parabolic boundaries:
INCLUDE_TOP_PARABOLA           = False
TOP_PARABOLA_TIP               = 1000e-9
TOP_PARABOLA_FOCUS             = 100e-9

INCLUDE_BOTTOM_PARABOLA        = True
BOTTOM_PARABOLA_TIP            = 0.0
BOTTOM_PARABOLA_FOCUS          = 300e-9
