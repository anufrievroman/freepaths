"""Config file to simulate a fishbone nanowire in Si at 4K"""

# General parameters:
OUTPUT_FOLDER_NAME             = 'Fishbone nanowire'
NUMBER_OF_PHONONS              = 3000
NUMBER_OF_TIMESTEPS            = 30000
T                              = 4.0


# Multiprocessing
NUMBER_OF_PROCESSES = 10


# Simulation time parameters:
TIMESTEP                       = 1.0e-12
total_simulation_time = 300e-9 # This should be at least a couple times the initialization time
NUMBER_OF_VIRTUAL_TIMESTEPS    = int(total_simulation_time / TIMESTEP)
initialization_time = 20e-9 # This should be set so that it is bigger than most phonons travel times
INITIALIZATION_TIMESTEPS       = int(initialization_time / TIMESTEP)
NUMBER_OF_INITIALIZATION_TIMEFRAMES = 3


# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 500e-9
LENGTH                         = 1500e-9


# Map & profiles parameters:
pixel_size = 10e-9
NUMBER_OF_PIXELS_X             = int(WIDTH / pixel_size)
NUMBER_OF_PIXELS_Y             = int(LENGTH / pixel_size)
IGNORE_FAULTY_PHONONS          = False


# Material parameters:
MEDIA                          = "Si"


# Internal scattering:
INCLUDE_INTERNAL_SCATTERING    = True
USE_GRAY_APPROXIMATION_MFP     = False
GRAY_APPROXIMATION_MFP         = 200e-9


# Phonon source:
PHONON_SOURCES                 = [Source(x=0, y=0, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="random_up")]


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
