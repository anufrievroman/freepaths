"""Test file which contains all possible shapes for testing"""


# General parameters:
OUTPUT_FOLDER_NAME             = 'All shapes'
NUMBER_OF_PHONONS              = 1000
NUMBER_OF_TIMESTEPS            = 30000
NUMBER_OF_NODES                = 400
TIMESTEP                       = 1.0e-12
T                              = 4.0
OUTPUT_TRAJECTORIES_OF_FIRST   = 30
NUMBER_OF_LENGTH_SEGMENTS      = 10


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
WIDTH                          = 700e-9
LENGTH                         = 2000e-9

# Phonon source:
PHONON_SOURCES                 = [Source(x=0, y=0, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="random_up")]

# Holes
HOLES                          = [
    CircularHole(x=-200e-9, y=200e-9, diameter=100e-9),
    RectangularHole(x=0, y=200e-9, size_x=100e-9, size_y=150e-9),
    TriangularUpHole(x=200e-9, y=200e-9, size_x=100e-9, size_y=150e-9),
    TriangularDownHole(x=-200e-9, y=400e-9, size_x=100e-9, size_y=150e-9),
    TriangularDownHalfHole(x=0, y=400e-9, size_x=100e-9, size_y=150e-9, is_right_half=False),
    TriangularDownHalfHole(x=200e-9, y=400e-9, size_x=100e-9, size_y=150e-9, is_right_half=True),
    TriangularUpHalfHole(x=-200e-9, y=600e-9, size_x=100e-9, size_y=150e-9, is_right_half=False),
    TriangularUpHalfHole(x=0, y=600e-9, size_x=100e-9, size_y=150e-9, is_right_half=True),
    ParabolaTop(tip=1000e-9, focus=100e-9),
]

# Roughness [m]:
SIDE_WALL_ROUGHNESS            = 2e-9
TOP_ROUGHNESS                  = 0.2e-9
BOTTOM_ROUGHNESS               = 0.2e-9

# Multiprocessing
NUMBER_OF_PROCESSES = 8
