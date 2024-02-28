"""Config file to simulate a structure with parabolic surface at the bottom
which acts like a collimator for phonon flux emitted from the hot spot of 10x10nm"""

# General parameters:
OUTPUT_FOLDER_NAME             = 'Parabolic lens collimation'
NUMBER_OF_PHONONS              = 2000
NUMBER_OF_TIMESTEPS            = 30000
T                              = 4.0


# Multiprocessing
NUMBER_OF_PROCESSES = 10

# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 2000e-9
LENGTH                         = 2000e-9


# Map & profiles parameters:
pixel_size = 30e-9
NUMBER_OF_PIXELS_X             = int(WIDTH / pixel_size)
NUMBER_OF_PIXELS_Y             = int(LENGTH / pixel_size)
IGNORE_FAULTY_PHONONS          = False


# Material parameters:
MEDIA                          = 'Si'


# Internal scattering:
INCLUDE_INTERNAL_SCATTERING    = True
USE_GRAY_APPROXIMATION_MFP     = False
GRAY_APPROXIMATION_MFP         = None


# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 2000e-9
LENGTH                         = 2000e-9


# Roughness [m]:
SIDE_WALL_ROUGHNESS            = 2e-9
TOP_ROUGHNESS                  = 0.2e-9
BOTTOM_ROUGHNESS               = 0.2e-9


# Hot and cold sides:
COLD_SIDE_POSITION_TOP         = True
COLD_SIDE_POSITION_LEFT        = True
COLD_SIDE_POSITION_RIGHT       = True

INCLUDE_RIGHT_SIDEWALL         = False
INCLUDE_LEFT_SIDEWALL          = False


# Phonon source:
PHONON_SOURCES = [Source(x=0, y=300e-9, z=0, size_x=100e-9,  size_y=100e-9, size_z=THICKNESS, angle_distribution="uniform")]


# Parabolic mirror:
HOLES = [ParabolaBottom(tip=0, focus=300e-9)]

