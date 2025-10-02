"""Config file to simulate a membrane with an array of interfaces"""

# General parameters:
OUTPUT_FOLDER_NAME             = 'Interfaces'
NUMBER_OF_PARTICLES            = 3000
T                              = 300
OUTPUT_TRAJECTORIES_OF_FIRST   = 10
LOW_MEMORY_USAGE               = True

# Simulation time parameters:
TIMESTEP                       = 0.5e-12
NUMBER_OF_TIMESTEPS            = 700000
NUMBER_OF_VIRTUAL_TIMESTEPS    = NUMBER_OF_TIMESTEPS*8
NUMBER_OF_STABILIZATION_TIMEFRAMES = 5
NUMBER_OF_TIMEFRAMES            = 8


# Multiprocessing:
NUMBER_OF_PROCESSES = 8

# Material parameters:
MEDIA                          = 'Si'

# Internal scattering:
INCLUDE_INTERNAL_SCATTERING    = True
USE_GRAY_APPROXIMATION_MFP     = False
GRAY_APPROXIMATION_MFP         = None

RETHERMALIZATION_ON_HOT_SIDES = True

# System dimensions [m]:
WIDTH                          = 180.5e-9 # security of 0.5nm to avoid numerical issues
THICKNESS                      = 1000e-9
LENGTH                         = 2000e-9



# Map & profiles parameters:
pixel_size =15e-9
NUMBER_OF_PIXELS_X             = int(WIDTH / pixel_size)
NUMBER_OF_PIXELS_Y             = int(LENGTH / pixel_size)
IGNORE_FAULTY_PARTICLES          = False

# PARTICLE source at the bottom:
PARTICLE_SOURCES                 = [Source(x=0, y=0, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="random", angle=0)]


# Interfaces (or thin layers of another material):
INTERFACES = []
period = 20e-9 # i.e. istance between interfaces
start_x = -WIDTH / 2 + period
end_x = WIDTH / 2 - period
INTERFACE_ROUGHNESS = 1.5e-9

x = start_x
while x <= end_x:
    INTERFACES.append(VerticalPlane(position_x=x, inner_material='Ge', outer_material=MEDIA, depth = THICKNESS))
    x += period




