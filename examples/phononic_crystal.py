"""Config file to simulate a phononic crystal with square lattice of holes"""

# General parameters:
OUTPUT_FOLDER_NAME             = 'Phononic crystal square'
NUMBER_OF_PHONONS              = 2000
NUMBER_OF_TIMESTEPS            = 30000
T                              = 4.0


# Multiprocessing:
NUMBER_OF_PROCESSES = 10


# Material parameters:
MEDIA                          = 'Si'


# Internal scattering:
INCLUDE_INTERNAL_SCATTERING    = True
USE_GRAY_APPROXIMATION_MFP     = False
GRAY_APPROXIMATION_MFP         = None


# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 1200e-9
LENGTH                         = 2200e-9


# Map & profiles parameters:
pixel_size = 15e-9
NUMBER_OF_PIXELS_X             = int(WIDTH / pixel_size)
NUMBER_OF_PIXELS_Y             = int(LENGTH / pixel_size)
IGNORE_FAULTY_PHONONS          = False


# Simulation time parameters:
TIMESTEP                       = 1.0e-12
total_simulation_time = 300e-9 # This should be at least a couple times the initialization time
NUMBER_OF_VIRTUAL_TIMESTEPS    = int(total_simulation_time / TIMESTEP)
initialization_time = 20e-9 # This should be set so that it is bigger than most phonons travel times
INITIALIZATION_TIMESTEPS       = int(initialization_time / TIMESTEP)
NUMBER_OF_INITIALIZATION_TIMEFRAMES = 3


# Phonon source:
PHONON_SOURCES                 = [Source(x=0, y=0, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="random_up")]

# Lattice of holes:
HOLES = []
period = 300e-9
for row in range(5):
    for column in range(5):
        x = - 4 * period / 2 + column * period
        y = (row + 1) * period
        HOLES.append(RectangularHole(x=x, y=y, size_x=200e-9, size_y=200e-9))
