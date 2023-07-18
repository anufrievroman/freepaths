"""Config file to simulate a structure with a parabolic surface at the top
which acts like a focusing mirror for a flux of parallel phonons emitted vertically"""

# General parameters:
OUTPUT_FOLDER_NAME             = 'Parabolic lens focusing'
NUMBER_OF_PHONONS              = 200
NUMBER_OF_TIMESTEPS            = 1000
NUMBER_OF_NODES                = 400
TIMESTEP                       = 1.0e-12
T                              = 4.0
PLOTS_IN_TERMINAL              = False
OUTPUT_SCATTERING_MAP          = False
OUTPUT_RAW_THERMAL_MAP         = True
OUTPUT_TRAJECTORIES_OF_FIRST   = 30
NUMBER_OF_LENGTH_SEGMENTS      = 10


# Map & profiles parameters:
NUMBER_OF_PIXELS_X             = 50
NUMBER_OF_PIXELS_Y             = 100
NUMBER_OF_TIMEFRAMES           = 6


# Material parameters:
MEDIA                          = 'Si'
SPECIFIC_HEAT_CAPACITY         = 0.0176  # [J/kg/K] for Si at 4 K


# Internal scattering:
INCLUDE_INTERNAL_SCATTERING    = True
USE_GRAY_APPROXIMATION_MFP     = False
GRAY_APPROXIMATION_MFP         = None


# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 1000e-9
LENGTH                         = 1100e-9


# Roughness [m]:
SIDE_WALL_ROUGHNESS            = 2e-9
TOP_ROUGHNESS                  = 0.2e-9
BOTTOM_ROUGHNESS               = 0.2e-9


# Phonon source:
PHONON_SOURCES                 = [Source(x=0, y=0, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="directional")]


# Parabolic boundary:
HOLES = [ParabolaTop(tip=1000e-9, focus=100e-9)]
