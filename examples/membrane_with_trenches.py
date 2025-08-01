"""Config file to simulate a membrane with array of trenches"""

# General parameters:
OUTPUT_FOLDER_NAME             = 'Membrane with trenches'
NUMBER_OF_PARTICLES            = 30
T                              = 300

# Simulation time parameters:
TIMESTEP                       = 0.5e-13
NUMBER_OF_TIMESTEPS            = 300000
NUMBER_OF_VIRTUAL_TIMESTEPS    = NUMBER_OF_TIMESTEPS*4
NUMBER_OF_STABILIZATION_TIMEFRAMES = 5
NUMBER_OF_TIMEFRAMES = 8

# Electron parameters: [eV]
ENERGY_UPPER_BOUND               = 100e-3
ENERGY_STEP                      = 15e-3
ELECTRON_MFP                     = 15e-9

# Multiprocessing:
NUMBER_OF_PROCESSES = 15

# Material parameters:
MEDIA                          = 'Si'

# Internal scattering:
INCLUDE_INTERNAL_SCATTERING    = True
USE_GRAY_APPROXIMATION_MFP     = False
GRAY_APPROXIMATION_MFP         = None

# System dimensions [m]:
THICKNESS                      = 300e-9
WIDTH                          = 3000e-9
LENGTH                         = 3000e-9

# Map & profiles parameters:
pixel_size = 15e-9
NUMBER_OF_PIXELS_X             = int(WIDTH / pixel_size)
NUMBER_OF_PIXELS_Y             = int(LENGTH / pixel_size)
IGNORE_FAULTY_PARTICLES        = False

# Particle source at the bottom:
PARTICLE_SOURCES               = [Source(x=0, y=0, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="random", angle=0)]

# Lattice of trenches (i.e. holes with finite depth):
HOLES = []
period = 500e-9
for row in range(5):

    # Perpendicular trenches:
    HOLES.append(RectangularHole(x=0, y=period+row*period, size_x=WIDTH, size_y=250e-9, depth=150e-9))

    # Parallel trenches:
    # HOLES.append(RectangularHole(x=-WIDTH/2+period+row*period, y=LENGTH/2, size_x=250e-9, size_y=LENGTH, depth=150e-9))
