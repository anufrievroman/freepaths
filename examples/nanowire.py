"""Config file to simulate a Si nanowire at 300 K"""


# General parameters:
OUTPUT_FOLDER_NAME             = 'Si nanowire at 300 K'
NUMBER_OF_PHONONS              = 10000
T                              = 300.0

# Simulation time parameters:
TIMESTEP                       = 1.0e-12
NUMBER_OF_TIMESTEPS            = 350000
NUMBER_OF_VIRTUAL_TIMESTEPS    = NUMBER_OF_TIMESTEPS*4
NUMBER_OF_STABILIZATION_TIMEFRAMES = 5
NUMBER_OF_TIMEFRAMES = 8

# Multiprocessing
NUMBER_OF_PROCESSES = 15

# Material parameters:
MEDIA                          = 'Si'

# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 200e-9
LENGTH                         = 2000e-9

# Map & profiles parameters:
pixel_size = 30e-9
NUMBER_OF_PIXELS_X             = int(WIDTH / pixel_size)
NUMBER_OF_PIXELS_Y             = int(LENGTH / pixel_size)
IGNORE_FAULTY_PHONONS          = False

# Phonon source:
PHONON_SOURCES                 = [Source(x=0, y=0, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="random", angle=0)]
