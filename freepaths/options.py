"""Module provides enums for the parameters"""

import enum


class Shapes(enum.Enum):
    """Possible shapes of holes"""
    CIRCLE = 1
    RECTANGLE = 2
    TRIANGLE_UP = 3
    TRIANGLE_DOWN = 4


class Materials(enum.Enum):
    """Possible materials"""
    SILICON = 1
    SILICON_CARBIDE = 2
    DIAMOND = 3
    ALUMINIUM_NITRIDE = 4


class Distributions(enum.Enum):
    """Possible distributions of angles"""
    RANDOM = 1
    LAMBERT = 2
    DIRECTIONAL = 3


class Positions(enum.Enum):
    """Possible positions of cold side"""
    TOP = 1
    RIGHT = 2
    TOP_AND_RIGHT = 3
    TOP_AND_BOTTOM = 4
