"""Module provides enums for the parameters"""

import enum


class Materials(enum.Enum):
    """Possible materials"""
    Si = 1
    SiC = 2
    Diamond = 3
    AlN = 4


class Distributions(enum.Enum):
    """Possible distributions of angles"""
    RANDOM = 1
    LAMBERT = 2
    DIRECTIONAL = 3
    UNIFORM = 4


class Positions(enum.Enum):
    """Possible positions of cold side"""
    TOP = 1
    RIGHT = 2
    TOP_AND_RIGHT = 3
    TOP_AND_BOTTOM = 4
