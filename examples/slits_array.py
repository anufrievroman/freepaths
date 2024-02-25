"""Config file to simulate a membrane with a staggered lattice of rectangular slits"""

import numpy as np


# General parameters:
OUTPUT_FOLDER_NAME             = 'Array of slits'
NUMBER_OF_PHONONS              = 2000
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


# Phonon source:
PHONON_SOURCES                 = [Source(x=0, y=0, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="random", angle=0)]


# Map & profiles parameters:
pixel_size = 15e-9
NUMBER_OF_PIXELS_X             = int(WIDTH / pixel_size)
NUMBER_OF_PIXELS_Y             = int(LENGTH / pixel_size)
IGNORE_FAULTY_PHONONS          = False


# Simulation time parameters:
TIMESTEP                       = 0.5e-12
total_simulation_time = 300e-9 # This should be at least a couple times the initialization time
NUMBER_OF_VIRTUAL_TIMESTEPS    = int(total_simulation_time / TIMESTEP)
initialization_time = 50e-9 # This should be set so that it is bigger than most phonons travel times
INITIALIZATION_TIMESTEPS       = int(initialization_time / TIMESTEP)
NUMBER_OF_INITIALIZATION_TIMEFRAMES = 3


# Staggered lattice of holes:
HOLES = []

size_x                         = 200e-9
size_y                         = 100e-9
period_x                       = 300e-9
period_y                       = 300e-9
first_hole_coordinate = 300e-9
number_of_periods_x = 5
number_of_periods_y = 3
for i in range(number_of_periods_y):
    for j in range(number_of_periods_x):
        x = -(number_of_periods_x - 1) * period_x / 2 + j * period_x
        y = first_hole_coordinate + i * period_y * 2
        HOLES.append(RectangularHole(x, y, size_x, size_y))
for i in range(number_of_periods_y):
    for j in range(number_of_periods_x - 1):
        x = - 2 * period_x / 2 + j * period_x - period_x / 2
        y = first_hole_coordinate + period_y + i * period_y * 2
        HOLES.append(RectangularHole(x, y, size_x, size_y))
