"""Test file which contains all possible shapes for testing"""


# General parameters:
OUTPUT_FOLDER_NAME             = 'All shapes'
NUMBER_OF_PHONONS              = 1000
NUMBER_OF_TIMESTEPS            = 30000
T                              = 4.0


# Multiprocessing
NUMBER_OF_PROCESSES = 10


# Simulation time parameters:
TIMESTEP                       = 1.0e-12
total_simulation_time = 300e-9 # This should be at least a couple times the initialization time
NUMBER_OF_VIRTUAL_TIMESTEPS    = int(total_simulation_time / TIMESTEP)
initialization_time = 50e-9 # This should be set so that it is bigger than most phonons travel times
INITIALIZATION_TIMESTEPS       = int(initialization_time / TIMESTEP)
NUMBER_OF_INITIALIZATION_TIMEFRAMES = 3


# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 700e-9
LENGTH                         = 1500e-9


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


# Phonon source:
PHONON_SOURCES                 = [Source(x=-WIDTH/2, y=LENGTH/2, z=0, size_x=0,  size_y=500e-9, size_z=THICKNESS, angle_distribution="random", angle=np.pi/2)]

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
