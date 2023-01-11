"""Module that provides phonon scattering types that occur on each step"""

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

    @property
    def is_diffuse(self):
        """Is any of the scattering types diffuse?"""
        return any([self.holes == Scattering.DIFFUSE,
                    self.pillars == Scattering.DIFFUSE,
                    self.top_bottom == Scattering.DIFFUSE,
                    self.walls == Scattering.DIFFUSE,
                    self.hot_side == Scattering.DIFFUSE])

    @property
    def is_internal(self):
        """Is any of the scattering types diffuse?"""
        return self.internal is not None

    @property
    def is_scattered(self):
        """Has any of the scattering events occurred?"""
        return any([self.holes is not None,
                    self.pillars is not None,
                    self.top_bottom is not None,
                    self.walls is not None,
                    self.internal is not None,
                    self.hot_side is not None])

    def reset(self):
        """Reset all scattering types to None"""
        self.holes = None
        self.pillars = None
        self.top_bottom = None
        self.walls = None
        self.internal = None
        self.hot_side = None
