"""Input file to simulate a 2D graphene sheet at 300 K"""

OUTPUT_FOLDER_NAME             = "Graphene sheet"
NUMBER_OF_PHONONS              = 50000
T                              = 300

# Simulation time parameters:
TIMESTEP                       = 1.0e-12
NUMBER_OF_TIMESTEPS            = 60000
NUMBER_OF_VIRTUAL_TIMESTEPS    = NUMBER_OF_TIMESTEPS*4
NUMBER_OF_STABILIZATION_TIMEFRAMES = 5
NUMBER_OF_TIMEFRAMES = 8

# Multiprocessing
NUMBER_OF_PROCESSES = 10

# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 1000e-9
LENGTH                         = 2200e-9

# Map & profiles parameters:
pixel_size = 50e-9
NUMBER_OF_PIXELS_X             = int(WIDTH / pixel_size)
NUMBER_OF_PIXELS_Y             = int(LENGTH / pixel_size)
IGNORE_FAULTY_PHONONS          = False

# Material:
MEDIA                          = "Graphite"
IS_TWO_DIMENSIONAL_MATERIAL    = True

PHONON_SOURCES = [Source(size_x=WIDTH, size_y=0, size_z=THICKNESS, angle_distribution="random", angle=0)]
