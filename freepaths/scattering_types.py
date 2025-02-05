"""
Module that provides types of phonon scattering that might occur.
These scattering types are returned on each step
so that higher level modules know what has happened.
"""

import enum


class Scattering(enum.Enum):
    """Possible scattering types"""
    DIFFUSE = 1
    SPECULAR = 2


class ScatteringTypes:
    """Phonon scattering types"""

    def __init__(self):
        """Initialize possible scattering type"""
        self.holes = None
        self.pillars = None
        self.top_bottom = None
        self.walls = None
        self.internal = None
        self.hot_side = None
        self.interfaces = None

    @property
    def is_diffuse(self):
        """Is any of the scattering types diffuse?"""
        return any([self.holes == Scattering.DIFFUSE,
                    self.pillars == Scattering.DIFFUSE,
                    self.top_bottom == Scattering.DIFFUSE,
                    self.walls == Scattering.DIFFUSE,
                    self.hot_side == Scattering.DIFFUSE,
                    self.interfaces == Scattering.DIFFUSE])

    @property
    def is_internal(self):
        """Is the scattering type internal?"""
        return self.internal is not None

    @property
    def is_diffuse_on_hole(self):
        """Was there a diffuse scattering on holes?"""
        return self.holes == Scattering.DIFFUSE

    @property
    def is_specular_on_hole(self):
        """Was there a specular scattering on holes?"""
        return self.holes == Scattering.SPECULAR

    @property
    def is_scattered(self):
        """Has any of the scattering events occurred?"""
        return any([self.holes,
                    self.pillars,
                    self.top_bottom,
                    self.walls,
                    self.internal,
                    self.hot_side,
                    self.interfaces])

    def reset(self):
        """Reset all scattering types to None"""
        self.holes = None
        self.pillars = None
        self.top_bottom = None
        self.walls = None
        self.internal = None
        self.hot_side = None
        self.interfaces = None


class ScatteringPlaces:
    """Phonon scattering places on a triangle"""

    def __init__(self):
        """Initialize possible scattering places"""
        self.right_wall = None
        self.left_wall = None
        self.floor = None

    @property
    def is_scattered(self):
        """Has any of the scattering events occurred?"""
        return any([self.right_wall,
                    self.left_wall,
                    self.floor])

    def reset(self):
        """Reset all scattering places to None"""
        self.right_wall = None
        self.left_wall = None
        self.floor = None
