"""
Module that provides various types of scattering objects like hole of different shapes.
User is expected to use one of these objects in the config file.
"""

import numpy


class CircularPillar:
    """Shape of a circular pillar with inclined wall"""

    def __init__(
        self, x=0, y=0, diameter=200e-9, height=300e-9, wall_angle=numpy.pi / 2
    ):
        self.x0 = x
        self.y0 = y
        self.diameter = diameter
        self.height = height
        self.wall_angle = wall_angle
