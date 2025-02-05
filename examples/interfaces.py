"""Config file to simulate a membrane with an array of interfaces"""

# General parameters:
OUTPUT_FOLDER_NAME             = 'Membrane with interfaces'
NUMBER_OF_PHONONS              = 1000
T                              = 300

# Simulation time parameters:
TIMESTEP                       = 0.5e-12
NUMBER_OF_TIMESTEPS            = 500000
NUMBER_OF_VIRTUAL_TIMESTEPS    = NUMBER_OF_TIMESTEPS*4
NUMBER_OF_STABILIZATION_TIMEFRAMES = 5
NUMBER_OF_TIMEFRAMES            = 8

# Multiprocessing:
NUMBER_OF_PROCESSES = 6

# Material parameters:
MEDIA                          = 'Si'

# Internal scattering:
INCLUDE_INTERNAL_SCATTERING    = True
USE_GRAY_APPROXIMATION_MFP     = False
GRAY_APPROXIMATION_MFP         = None

# System dimensions [m]:
THICKNESS                      = 1500e-9
WIDTH                          = 3000e-9
LENGTH                         = 3000e-9

# Map & profiles parameters:
pixel_size = 15e-9
NUMBER_OF_PIXELS_X             = int(WIDTH / pixel_size)
NUMBER_OF_PIXELS_Y             = int(LENGTH / pixel_size)
IGNORE_FAULTY_PHONONS          = False

# Phonon source at the bottom:
PHONON_SOURCES                 = [Source(x=0, y=0, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="random", angle=0)]

# Lattice of interfaces
INTERFACES = []
period = 300e-9
for index in range(9):
    INTERFACES.append(VerticalPlane(position_x=-4*period+period*index, transmission=0.5))
    # INTERFACES.append(HorizontalPlane(position_z=-2*period+period*index, transmission=0.5))

