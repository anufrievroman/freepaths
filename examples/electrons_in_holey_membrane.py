"""Config file to simulate a membrane with array of holes"""

# General parameters:
OUTPUT_FOLDER_NAME             = 'Holey membrane'
NUMBER_OF_PARTICLES            = 1000 * 8
T                              = 300

# Simulation time parameters:
TIMESTEP                       = 1e-14
NUMBER_OF_TIMESTEPS            = 300000
NUMBER_OF_VIRTUAL_TIMESTEPS    = NUMBER_OF_TIMESTEPS*4
NUMBER_OF_STABILIZATION_TIMEFRAMES = 6
NUMBER_OF_TIMEFRAMES = 8

# Electron parameters: [eV]
ENERGY_UPPER_BOUND               = 160e-3
ENERGY_STEP                      = 10e-3
ELECTRON_MFP                     = 15e-9
MEAN_MAPPING_CONSTANT            = 1.9e-5

# Multiprocessing:
NUMBER_OF_PROCESSES = 15

# Material parameters:
from scipy.constants import electron_volt
MEDIA                          = 'Si'
MEDIA_FERMI_LEVEL              = -0.3* electron_volt

# Internal scattering:
INCLUDE_INTERNAL_SCATTERING    = True
USE_GRAY_APPROXIMATION_MFP     = False
GRAY_APPROXIMATION_MFP         = None

# System dimensions [m]:
period = 60e-9
neck = 30e-9
THICKNESS                      = 1e-6
WIDTH                          = period * 4
LENGTH                         = period * 4

# Map & profiles parameters:
pixel_size = 3e-9
NUMBER_OF_PIXELS_X             = int(WIDTH / pixel_size)
NUMBER_OF_PIXELS_Y             = int(LENGTH / pixel_size)
IGNORE_FAULTY_PARTICLES        = False

# Particle source at the bottom:
PARTICLE_SOURCES               = [Source(x=0, y=0, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="random", angle=0)]

# Hexagonal lattice of circular holes:
HOLES = []
rows, cols = 8,5
for row in range(rows):
    for col in range(cols):
        x = (col+(1/2)*(row%2==0)) * period - (cols-cols//2 - 1/2) * period
        y = (row+1) * period / 2
        if row%2==1 and col==0:
            continue
        diam = ((2**0.5)/2) * period - neck
        HOLES.append(CircularHole(x=x, y=y, diameter=diam))
