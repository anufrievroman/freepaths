"""Example file to explain the point line hole shape"""


# General parameters:
OUTPUT_FOLDER_NAME             = 'Point line example'
NUMBER_OF_PHONONS              = 1000
NUMBER_OF_TIMESTEPS            = 30000
T                              = 4.0


# Multiprocessing
NUMBER_OF_PROCESSES = 10


# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 700e-9
LENGTH                         = 1500e-9


# Map & profiles parameters:
pixel_size = 10e-9
NUMBER_OF_PIXELS_X             = int(WIDTH / pixel_size)
NUMBER_OF_PIXELS_Y             = int(LENGTH / pixel_size)
IGNORE_FAULTY_PHONONS          = False


# Sources
PHONON_SOURCES                 = [Source(x=0, y=0, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="random_right")]


# Holes
points                         = [(0e-9, 0e-9), (50e-9, 100e-9)]
HOLES                          = [
    PointLineHole(x=0e-9, y=0e-9, points=points),
    PointLineHole(x=50e-9, y=200e-9, points=points),
    ]
