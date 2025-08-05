"""This module provides electron class which generates and moves an electron"""

from random import choice
from scipy.constants import h, electron_volt
from freepaths.config import cf
from freepaths.particle import Particle
from freepaths.particle_types import ParticleType

import numpy as np

class Electron(Particle):
    """An electron particle"""
    
    def __init__(self, material):
        super().__init__()
        
        # Assign particle type
        self.type = ParticleType.ELECTRON 
        self.is_electron_carrier = cf.is_carrier_electron # assign electron or hole behavior
        
        source = choice(cf.particles_sources)
        while True:
            self.x, self.y, self.z = source.generate_coordinates()
            is_in_hole = any(hole.is_inside(self.x, self.y, None, cf) for hole in cf.holes)
            if not is_in_hole:
                break
        
        # Assign initial angles
        self.theta, self.phi = source.generate_angles()
        self.correct_angle()
        if cf.is_two_dimensional_material:
            self.phi, self.z = 0.0, 0.0
            
        self.assign_energy()
        self.assign_frequency(material)
        self.assign_speed(material)
        self.assign_internal_scattering_time(material)
    
    def assign_energy(self):
        """
        Assign energy uniformly to an electron based on temperature.
        Conduction minimum is considered at 0 energy by convention.
        """
        num_energy_points = int((cf.energy_upper_bound-cf.energy_lower_bound)/cf.energy_step)
        # Keep energy in J for other computations
        self.energy = np.random.choice(np.linspace(cf.energy_lower_bound, cf.energy_upper_bound - cf.energy_step, num_energy_points)) * electron_volt + cf.energy_step * electron_volt
        
    
    @property
    def wavelength(self):
        raise Exception("Electron wavelength undefined")
    
    def assign_frequency(self, material):
        """
        Conduction minimum is 0 by convention
        """
        self.f = self.energy / h
    
    def assign_speed(self, material):
        """
        Assign speed to an electron using group velocity and considering parabolic energy-k relation and conduction minimum at 0
        """
        if self.is_electron_carrier:
            effective_mass = material.effective_electron_dos_mass
        else:
            effective_mass = material.effective_hole_dos_mass
        self.speed = (2*self.energy/(effective_mass))**(0.5)
    
    def assign_internal_scattering_time(self, material):
        """
        Determine relaxation time after which this particle will undergo internal scattering.
        Considering only elasctic scattering between electrons and phonons and thus a constant MFP.
        """
        mfp = cf.electron_mfp
        self.time_of_internal_scattering = mfp / self.speed
        
        # self.time_of_internal_scattering = -np.log(random()) * (mfp/self.speed)