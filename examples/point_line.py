"""Example file to explain the point line hole shape"""

import math
import numpy as np


# General parameters:
OUTPUT_FOLDER_NAME             = 'Point line example'
NUMBER_OF_PHONONS              = 100
NUMBER_OF_TIMESTEPS            = 30000
T                              = 100
TIMESTEP                       = 1e-12


# Multiprocessing
NUMBER_OF_PROCESSES            = 10
HOLE_ROUGHNESS                 = 0.2e-9


# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 700e-9
LENGTH                         = 700e-9


# Map & profiles parameters:
pixel_size = 5e-9
NUMBER_OF_PIXELS_X             = int(WIDTH / pixel_size)
NUMBER_OF_PIXELS_Y             = int(LENGTH / pixel_size)
IGNORE_FAULTY_PHONONS          = False


# Sources
PHONON_SOURCES                 = [Source(x=0, y=0, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="random_up")]


# Holes

# Write a python function that takes radius, start_angle, end_angle, resolution and returns a list of points which are on an arc with center 0,0 respective radius and angles. the points have the distance specified in resolution

def points_on_arc(radius, start_angle, end_angle, resolution):
    # Convert angles to radians
    start_angle_rad = math.radians(start_angle)
    end_angle_rad = math.radians(end_angle)

    # Calculate the number of points
    num_points = int(radius* abs(end_angle_rad - start_angle_rad) / resolution) + 1

    # Calculate angular increment
    angle_increment = (end_angle_rad - start_angle_rad) / (num_points - 1)

    # Generate points on the arc
    arc_points = []
    for i in range(num_points):
        angle = start_angle_rad + i * angle_increment
        x = radius * math.sin(angle)
        y = radius * math.cos(angle)
        arc_points.append((x, y))

    return arc_points

def make_arc_hole(radius, start_angle, end_angle, resolution, x, y, thickness):
    points = points_on_arc(radius, start_angle, end_angle, resolution)
    return PointLineHole(x=x, y=y, thickness=thickness, points=points)

sinxdivx = lambda x: 1 if x == 0 else np.sin(x)/x

HOLES = [
    PointLineHole(x=0, y=200e-9, points=[(0,0), (0e-9, 50e-9)], thickness=30e-9, rotation=0),
    PointLineHole(x=0, y=200e-9, points=[(0,0), (0e-9, 50e-9)], thickness=30e-9, rotation=45),
    ]

