"""Config file to simulate a membrane with a staggered lattice of rectangular slits"""

import numpy as np


# General parameters:
OUTPUT_FOLDER_NAME             = 'Array of triangles'
NUMBER_OF_PARTICLES            = 2000
TIMESTEP                       = 0.5e-12
NUMBER_OF_TIMESTEPS            = 30000
T                              = 4.0

# Multiprocessing
NUMBER_OF_PROCESSES = 10

# Material parameters:
MEDIA                          = 'Si'

# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 1200e-9
LENGTH                         = 2200e-9

# Particle source:
PARTICLE_SOURCES                 = [Source(x=0, y=0, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="random", angle=0)]

# Map & profiles parameters:
pixel_size = 15e-9
NUMBER_OF_PIXELS_X             = int(WIDTH / pixel_size)
NUMBER_OF_PIXELS_Y             = int(LENGTH / pixel_size)
IGNORE_FAULTY_PARTICLES        = False

# Staggered lattice of triangular holes:
HOLES = []

period_x = 350e-9
period_y = 350e-9
first_hole_coordinate = 300e-9
number_of_periods_x = 3
number_of_periods_y = 3

size_x = 100e-9
size_y = 200e-9
for i in range(number_of_periods_y):
    for j in range(number_of_periods_x):
        x = -(number_of_periods_x - 1) * period_x / 2 + j * period_x
        y = first_hole_coordinate + i * period_y * 2
        HOLES.append(TriangularUpHole(x, y, size_x*(0.4*i + 1), size_y*(0.4*i + 1)))
for i in range(number_of_periods_y):
    for j in range(number_of_periods_x + 1):
        x = - 2 * period_x / 2 + j * period_x - period_x / 2
        y = first_hole_coordinate + period_y + i * period_y * 2
        HOLES.append(TriangularDownHole(x, y, size_x*(0.4*i + 1), size_y*(0.4*i + 1)))
