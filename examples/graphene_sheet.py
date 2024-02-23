"""Simple input file to simulate a 2D graphene sheet"""

OUTPUT_FOLDER_NAME             = "Graphene sheet"
NUMBER_OF_PHONONS              = 10000
NUMBER_OF_TIMESTEPS            = 60000
T                              = 300


# Multiprocessing
NUMBER_OF_PROCESSES = 10


# Simulation time parameters:
TIMESTEP                       = 1e-12
total_simulation_time = 200e-9 # This should be at least a couple times the initialization time
NUMBER_OF_VIRTUAL_TIMESTEPS    = int(total_simulation_time / TIMESTEP)
initialization_time = 20e-9 # This should be set so that it is bigger than most phonons travel times
INITIALIZATION_TIMESTEPS       = int(initialization_time / TIMESTEP)
NUMBER_OF_INITIALIZATION_TIMEFRAMES = 3


# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 1000e-9
LENGTH                         = 2200e-9


# Map & profiles parameters:
pixel_size = 30e-9
NUMBER_OF_PIXELS_X             = int(WIDTH / pixel_size)
NUMBER_OF_PIXELS_Y             = int(LENGTH / pixel_size)
IGNORE_FAULTY_PHONONS          = False


# Material:
MEDIA                          = "Graphite"
IS_TWO_DIMENSIONAL_MATERIAL    = True
SPECIFIC_HEAT_CAPACITY         = 710.0  # [J/kg/K] for Graphene at 300 K

PHONON_SOURCES = [Source(size_x=WIDTH, size_y=0, size_z=THICKNESS, angle_distribution="random_up")]
