"""
Module that provides various types of scattering objects like hole of different shapes.
User is expected to use one of these objects in the config file.
"""

import numpy

class CircularHole:
    """Shape of a circular hole"""
    def __init__(self, x=0, y=0, diameter=100e-9):
        self.x = x
        self.y = y
        self.diameter = diameter


class RectangularHole:
    """Shape of a rectangular hole"""
    def __init__(self, x=0, y=0, size_x=100e-9, size_y=100e-9):
        self.x = x
        self.y = y
        self.size_x = size_x
        self.size_y = size_y


class TriangularUpHole:
    """Shape of a triangular hole facing up"""
    def __init__(self, x=0, y=0, size_x=100e-9, size_y=100e-9):
        self.x = x
        self.y = y
        self.size_x = size_x
        self.size_y = size_y


class TriangularDownHole:
    """Shape of a triangular hole facing down"""
    def __init__(self, x=0, y=0, size_x=100e-9, size_y=100e-9):
        self.x = x
        self.y = y
        self.size_x = size_x
        self.size_y = size_y


class TriangularDownHalfHole:
    """Shape of a half triangular hole facing down"""
    def __init__(self, x=0, y=0, size_x=100e-9, size_y=100e-9, is_right_half=True):
        self.x = x
        self.y = y
        self.size_x = size_x
        self.size_y = size_y
        self.is_right_half = is_right_half


class TriangularUpHalfHole:
    """Shape of a half triangular hole facing up"""
    def __init__(self, x=0, y=0, size_x=100e-9, size_y=100e-9, is_right_half=True):
        self.x = x
        self.y = y
        self.size_x = size_x
        self.size_y = size_y
        self.is_right_half = is_right_half


class CircularPillar:
    """Shape of a circular pillar with inclined wall"""
    def __init__(self, x=0, y=0, diameter=200e-9, height=300e-9, wall_angle=numpy.pi/2):
        self.x = x
        self.y = y
        self.diameter = diameter
        self.height = height
        self.wall_angle = wall_angle


class ParabolaTop:
    """Shape of a parabolic wall"""
    def __init__(self, tip=0, focus=0):
        self.tip = tip
        self.focus = focus


class ParabolaBottom:
    """Shape of a parabolic wall"""
    def __init__(self, tip=0, focus=0):
        self.tip = tip
        self.focus = focus
