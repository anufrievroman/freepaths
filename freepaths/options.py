"""Module provides enums for the parameters"""

import enum


class Materials(enum.Enum):
    """Possible materials"""
    Si = 1
    SiC = 2
    Diamond = 3
    AlN = 4
    Graphite = 5


class Distributions(enum.Enum):
    """Possible distributions of angles"""
    RANDOM_UP = 1
    RANDOM_DOWN = 2
    RANDOM_RIGHT = 3
    RANDOM_LEFT = 4
    LAMBERT = 5
    DIRECTIONAL = 6
    DIRECTIONAL_DOWN = 7
    UNIFORM = 8
