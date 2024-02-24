"""Config file to simulate a structure with a parabolic surface at the top
which acts like a focusing mirror for a flux of parallel phonons emitted vertically"""

# General parameters:
OUTPUT_FOLDER_NAME             = 'Parabolic lens focusing'
NUMBER_OF_PHONONS              = 1000
NUMBER_OF_TIMESTEPS            = 10000
T                              = 4.0


# Multiprocessing
NUMBER_OF_PROCESSES = 10


# Material parameters:
MEDIA                          = 'Si'


# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 1000e-9
LENGTH                         = 1100e-9


# Map & profiles parameters:
pixel_size = 10e-9
NUMBER_OF_PIXELS_X             = int(WIDTH / pixel_size)
NUMBER_OF_PIXELS_Y             = int(LENGTH / pixel_size)
IGNORE_FAULTY_PHONONS          = False


# Adjust walls
COLD_SIDE_POSITION_BOTTOM        = True
HOT_SIDE_POSITION_BOTTOM         = False


# Phonon source:
PHONON_SOURCES                 = [Source(x=0, y=10e-9, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="directional")]


# Parabolic boundary:
HOLES = [ParabolaTop(tip=1000e-9, focus=100e-9)]
