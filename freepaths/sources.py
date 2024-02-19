"""Module provides the phonon source object"""

from random import random
from math import pi, asin
from numpy import sign
from freepaths.options import Distributions


class Source:
    """Phonon source as rectangular are where phonons are initiated"""
    def __init__(self, x=0, y=0, z=0, size_x=0, size_y=0, size_z=0, angle_distribution="random_up"):
        self.x = x
        self.y = y
        self.z = z
        self.size_x = size_x
        self.size_y = size_y
        self.size_z = size_z
        self.angle_distribution = angle_distribution

    def generate_coordinates(self):
        """Generate coordinates of the phonon inside the source"""
        phonon_x = self.x + 0.49 * self.size_x * (2 * random() - 1)
        phonon_y = self.y + 0.49 * self.size_y * (2 * random() - 1)
        phonon_z = self.z + 0.49 * self.size_z * (2 * random() - 1)
        return phonon_x, phonon_y, phonon_z

    def generate_angles(self):
        """Generate angles of the phonon inside the source"""

        if self.angle_distribution == Distributions.RANDOM_UP:
            theta = -pi/2 + pi*random()
            phi = asin(2*random() - 1)
            return theta, phi

        if self.angle_distribution == Distributions.RANDOM_DOWN:
            rand_sign = sign((2*random() - 1))
            theta = rand_sign*(pi/2 + pi/2*random())
            phi = asin(2*random() - 1)
            return theta, phi

        if self.angle_distribution == Distributions.RANDOM_RIGHT:
            theta = pi*random()
            phi = asin(2*random() - 1)
            return theta, phi

        if self.angle_distribution == Distributions.RANDOM_LEFT:
            theta = - pi*random()
            phi = asin(2*random() - 1)
            return theta, phi

        if self.angle_distribution == Distributions.DIRECTIONAL:
            theta = 1e-10
            phi = -pi/2 + pi*random()
            return theta, phi

        if self.angle_distribution == Distributions.LAMBERT:
            theta = asin(2*random() - 1)
            phi = asin((asin(2*random() - 1))/(pi/2))
            return theta, phi

        if self.angle_distribution == Distributions.UNIFORM:
            theta = -pi + 2*pi*random()
            phi = asin(2*random() - 1)
            return theta, phi
