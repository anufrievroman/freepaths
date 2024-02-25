"""Module provides the phonon source object"""

import enum

from random import random
from math import pi, asin


class Distributions(enum.Enum):
    """Possible distributions of angles"""
    RANDOM = 1
    LAMBERT = 2
    DIRECTIONAL = 3
    UNIFORM = 4


class Source:
    """Phonon source as rectangular are where phonons are initiated"""
    def __init__(self, x=0, y=0, z=0, size_x=0, size_y=0, size_z=0, angle_distribution="random_up", angle=0):
        self.x = x
        self.y = y
        self.z = z
        self.size_x = size_x
        self.size_y = size_y
        self.size_z = size_z
        self.angle_distribution = angle_distribution
        self.angle = angle

    def generate_coordinates(self):
        """Generate coordinates of the phonon inside the source"""
        phonon_x = self.x + 0.49 * self.size_x * (2 * random() - 1)
        phonon_y = self.y + 0.49 * self.size_y * (2 * random() - 1)
        phonon_z = self.z + 0.49 * self.size_z * (2 * random() - 1)
        return phonon_x, phonon_y, phonon_z

    def generate_angles(self):
        """Generate angles of the phonon inside the source"""
        if self.angle_distribution == Distributions.RANDOM:
            theta = -pi/2 + pi*random() + self.angle
            phi = asin(2*random() - 1)
        elif self.angle_distribution == Distributions.DIRECTIONAL:
            theta = self.angle + 1e-10
            phi = -pi/2 + pi*random()
        elif self.angle_distribution == Distributions.LAMBERT:
            theta = asin(2*random() - 1) + self.angle
            phi = asin((asin(2*random() - 1))/(pi/2))
        elif self.angle_distribution == Distributions.UNIFORM:
            theta = -pi + 2*pi*random()
            phi = asin(2*random() - 1)
        else:
            raise ValueError("Invalid distribution type")

        return theta, phi
