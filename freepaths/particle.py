"""This module provides an abstract class for particles"""

from abc import ABC, abstractmethod
from numpy import sign
from math import pi
from random import randint

from freepaths.config import cf
import freepaths.move
import numpy as np


class Particle(ABC):
    """A general particle for a MC simulation"""
    
    def __init__(self):
        # Type of the particle
        self.type = None
        
        # Position, geometric informations
        self.x = None
        self.y = None
        self.z = None
        self.vx = None
        self.vy = None
        self.vz = None
        self.phi = None
        self.theta = None
        
        # Physical properties
        self.speed = None
        self.f = None
        self.time_of_internal_scattering = None
        
        # Assign random first timestep
        self.first_timestep = randint(0, cf.number_of_virtual_timesteps)
    
    def update_speed_component(self):
        """Projette la vitesse scalaire sur vx, vy, vz selon theta/phi."""
        # Composantes sphériques
        sin_t = np.sin(self.theta)
        self.vx = self.speed * sin_t * np.cos(self.phi)
        self.vy = self.speed * sin_t * np.sin(self.phi)
        # En 2D, phi=0 → vz = 0
        self.vz = self.speed * np.cos(self.theta) if not cf.is_two_dimensional_material else 0.0
    
    
    @property
    @abstractmethod
    def wavelength(self):
        """Calculate wavelength of the particle"""
        pass
    
    @property
    def has_crossed_cold_side(self):
        """
        Checks if the particle at this timestep crossed the cold side.
        Depending on where user set cold sides, we check if particle crossed that line.
        Return boolean of wheather any of the cold sides has been crossed.
        """
        has_crossed_top = self.y > cf.length
        has_crossed_bottom = self.y < 0
        has_crossed_right = self.x > cf.width / 2.0
        has_crossed_left = self.x < - cf.width / 2.0
        return ((cf.cold_side_position_top and has_crossed_top) or
                (cf.cold_side_position_bottom and has_crossed_bottom) or
                (cf.cold_side_position_right and has_crossed_right) or
                (cf.cold_side_position_left and has_crossed_left))
        
    @property
    def has_crossed_hot_side(self):
        """
        Checks if the particle at this timestep crossed the hot side.
        Depending on where user set hot sides, we check if particle crossed that line.
        Return boolean of wheather any of the hot sides has been crossed.
        """
        has_crossed_top = self.y > cf.length
        has_crossed_bottom = self.y < 0
        has_crossed_right = self.x > cf.width / 2.0
        has_crossed_left = self.x < - cf.width / 2.0
        return ((cf.hot_side_position_top and has_crossed_top) or
                (cf.hot_side_position_bottom and has_crossed_bottom) or
                (cf.hot_side_position_right and has_crossed_right) or
                (cf.hot_side_position_left and has_crossed_left))
    
    def correct_angle(self):
        """Check if angles are out of the [-pi:pi] range and return them back to this range"""
        if abs(self.theta) > pi:
            self.theta -= sign(self.theta)*2*pi
    
    @abstractmethod
    def assign_frequency(self, material):
        """Assigning a frequency to the particle"""
        pass
    
    @abstractmethod
    def assign_speed(self, material):
        """Assining speed to the particle"""
        pass
    
    @abstractmethod
    def assign_internal_scattering_time(self, material):
        """Determine relaxation time after which this particle will undergo internal scattering"""
        pass
    
    def move(self):
        """Move the particle in one timestep and return new coordinates"""
        self.x, self.y, self.z = freepaths.move.move(self, cf.timestep)
        self.update_speed_component()