"""Example file to explain the point line hole shape"""

import math
import numpy as np


# General parameters:
OUTPUT_FOLDER_NAME             = 'Point line example'
NUMBER_OF_PARTICLES            = 100
NUMBER_OF_TIMESTEPS            = 30000
T                              = 50
TIMESTEP                       = 1e-12

# Multiprocessing
NUMBER_OF_PROCESSES            = 10

# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 300e-9
LENGTH                         = 300e-9

# Map & profiles parameters:
pixel_size = 5e-9
NUMBER_OF_PIXELS_X             = int(WIDTH / pixel_size)
NUMBER_OF_PIXELS_Y             = int(LENGTH / pixel_size)
IGNORE_FAULTY_PARTICLES        = False

# Sources
PARTICLE_SOURCES               = [Source(x=0, y=0, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="random", angle=0)]

# Holes

# See the Creating new holes (the easy way) tutorial in the documentation for more examples
def points_on_bezier(control_point1, handle1, control_point2, handle2, number_of_points):
    # get coordinate points
    x1, y1 = control_point1
    x2 = x1 + handle1[0]
    y2 = y1 + handle1[1]
    x4, y4 = control_point2
    x3 = x4 + handle2[0]
    y3 = y4 + handle2[1]
    # define function
    bezier_function = lambda t, p1, p2, p3, p4: (1-t)**3*p1 + 3*t*(1-t)**2*p2 + 3*t**2*(1-t)*p3 + t**3*p4
    bezier_point_function = lambda t: (bezier_function(t, x1, x2, x3, x4), bezier_function(t, y1, y2, y3, y4))

    # calculate points and return
    return [bezier_point_function(t) for t in np.linspace(0, 1, number_of_points)]

points = points_on_bezier((-100e-9, 0), (50e-9, 100e-9), (100e-9, 0), (-50e-9, -100e-9), 100)
HOLES = [PointLineHole(0, 150e-9, points, thickness=30e-9)]
