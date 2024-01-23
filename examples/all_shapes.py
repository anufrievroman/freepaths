"""Test file which contains all possible shapes for testing"""


# General parameters:
OUTPUT_FOLDER_NAME             = 'All shapes c'
NUMBER_OF_PHONONS              = 1000
NUMBER_OF_TIMESTEPS            = 30000
NUMBER_OF_VIRTUAL_TIMESTEPS    = 100000
NUMBER_OF_NODES                = 400
TIMESTEP                       = 1.0e-12
T                              = 4.0
OUTPUT_TRAJECTORIES_OF_FIRST   = 100
NUMBER_OF_LENGTH_SEGMENTS      = 10
IGNORE_FAULTY_PHONONS          = False


# Map & profiles parameters:
NUMBER_OF_PIXELS_X             = 200
NUMBER_OF_PIXELS_Y             = 200
NUMBER_OF_TIMEFRAMES           = 6


# Material parameters:
MEDIA                          = 'Si'
SPECIFIC_HEAT_CAPACITY         = 0.0176  # [J/kg/K] for Si at 4 K
#SPECIFIC_HEAT_CAPACITY         = 714     # [J/kg/K] for Si at 300 K


# Internal scattering:
INCLUDE_INTERNAL_SCATTERING    = True
USE_GRAY_APPROXIMATION_MFP     = False
GRAY_APPROXIMATION_MFP         = None


# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 700e-9
LENGTH                         = 1500e-9

# Phonon source:
PHONON_SOURCES                 = [Source(x=-WIDTH/2, y=LENGTH/2, z=0, size_x=0,  size_y=500e-9, size_z=THICKNESS, angle_distribution="random_up")]

# Walls:
INCLUDE_RIGHT_SIDEWALL           = False
INCLUDE_LEFT_SIDEWALL            = False
INCLUDE_TOP_SIDEWALL             = False
INCLUDE_BOTTOM_SIDEWALL          = False

# Cold sides:
COLD_SIDE_POSITION_TOP           = False
COLD_SIDE_POSITION_BOTTOM        = False
COLD_SIDE_POSITION_RIGHT         = True
COLD_SIDE_POSITION_LEFT          = False

# Hot sides:
HOT_SIDE_POSITION_TOP            = False
HOT_SIDE_POSITION_BOTTOM         = False
HOT_SIDE_POSITION_RIGHT          = False
HOT_SIDE_POSITION_LEFT           = True

# Holes
HOLES                          = [
    CircularHole(x=-200e-9, y=400e-9, diameter=100e-9),
    RectangularHole(x=0, y=400e-9, size_x=100e-9, size_y=150e-9),
    TriangularUpHole(x=200e-9, y=400e-9, size_x=100e-9, size_y=150e-9),
    TriangularDownHole(x=-200e-9, y=600e-9, size_x=100e-9, size_y=150e-9),
    TriangularDownHalfHole(x=0, y=600e-9, size_x=100e-9, size_y=150e-9, is_right_half=False),
    TriangularDownHalfHole(x=200e-9, y=600e-9, size_x=100e-9, size_y=150e-9, is_right_half=True),
    TriangularUpHalfHole(x=-200e-9, y=800e-9, size_x=100e-9, size_y=150e-9, is_right_half=False),
    TriangularUpHalfHole(x=0, y=800e-9, size_x=100e-9, size_y=150e-9, is_right_half=True),
    ParabolaTop(tip=LENGTH-100e-9, focus=100e-9),
    ParabolaBottom(tip=100e-9, focus=100e-9),
]

# Roughness [m]:
SIDE_WALL_ROUGHNESS            = 2e-9
TOP_ROUGHNESS                  = 0.2e-9
BOTTOM_ROUGHNESS               = 0.2e-9

# Multiprocessing
NUMBER_OF_PROCESSES = 50
