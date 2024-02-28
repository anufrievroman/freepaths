"""Config file to simulate a membrane with square lattice of pillars on the surface"""

# General parameters:
OUTPUT_FOLDER_NAME             = 'Membrane with pillars'
NUMBER_OF_PHONONS              = 500
T                              = 4.0

# Simulation time parameters:
TIMESTEP                       = 1.0e-12
NUMBER_OF_TIMESTEPS            = 30000
NUMBER_OF_VIRTUAL_TIMESTEPS    = NUMBER_OF_TIMESTEPS*4
NUMBER_OF_STABILIZATION_TIMEFRAMES = 5
NUMBER_OF_TIMEFRAMES = 8

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

# Lattice of pillars:
PILLARS = []
period = 300e-9
for row in range(3):
    for column in range(3):
        x = - 2 * period / 2 + column * period
        y = (row + 1) * period
        PILLARS.append(CircularPillar(x=x, y=y, diameter=200e-9, height=200e-9))
