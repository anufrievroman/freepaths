"""Config file to simulate a membrane with square lattice of pillars on the surface"""

# General parameters:
OUTPUT_FOLDER_NAME             = 'Membrane with pillars'
NUMBER_OF_PHONONS              = 500
NUMBER_OF_TIMESTEPS            = 30000
T                              = 4.0


# Multiprocessing
NUMBER_OF_PROCESSES = 10


# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 900e-9
LENGTH                         = 1200e-9


# Map & profiles parameters:
pixel_size = 30e-9
NUMBER_OF_PIXELS_X             = int(WIDTH / pixel_size)
NUMBER_OF_PIXELS_Y             = int(LENGTH / pixel_size)
IGNORE_FAULTY_PHONONS          = False


# Simulation time parameters:
TIMESTEP                       = 1.0e-12
total_simulation_time = 200e-9 # This should be at least a couple times the initialization time
NUMBER_OF_VIRTUAL_TIMESTEPS    = int(total_simulation_time / TIMESTEP)
initialization_time = 20e-9 # This should be set so that it is bigger than most phonons travel times
INITIALIZATION_TIMESTEPS       = int(initialization_time / TIMESTEP)
NUMBER_OF_INITIALIZATION_TIMEFRAMES = 3


# Material parameters:
MEDIA                          = 'Si'


# Internal scattering:
INCLUDE_INTERNAL_SCATTERING    = True
USE_GRAY_APPROXIMATION_MFP     = False
GRAY_APPROXIMATION_MFP         = None


# Phonon source:
PHONON_SOURCES                 = [Source(x=0, y=0, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="random", angle=0)]


# Roughness [m]:
SIDE_WALL_ROUGHNESS            = 2e-9
HOLE_ROUGHNESS                 = 2e-9
TOP_ROUGHNESS                  = 0.2e-9
BOTTOM_ROUGHNESS               = 0.2e-9


# Lattice of holes:
PILLARS = []
period = 300e-9
for row in range(3):
    for column in range(3):
        x = - 2 * period / 2 + column * period
        y = (row + 1) * period
        PILLARS.append(CircularPillar(x=x, y=y, diameter=200e-9, height=200e-9))
