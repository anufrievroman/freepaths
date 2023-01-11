"""This module provides phonon class which generates and moves a phonon"""

from math import pi, asin, exp, log
from random import random, choice
from numpy import sign
from scipy.constants import k, hbar
import numpy as np

import parameters as par
import move

class Phonon:
    """A phonon particle with various physical properties"""

    def __init__(self, material):
        """Initialize a phonon by assigning coordinates and other properties"""
        self.x = None
        self.y = None
        self.z = None
        self.f = None
        self.phi = None
        self.theta = None
        self.speed = None
        self.polarization = choice(['TA', 'TA', 'LA'])
        self.assign_initial_coordinates()
        self.assign_angles()
        self.assign_frequency(material)
        self.assign_speed(material)
        self.assign_internal_scattering_time(material)

    @property
    def wavelength(self):
        """Calculate wavelength on the phonon"""
        return self.speed/self.f

    @property
    def is_in_system(self):
        """Checks if the phonon at this timestep did not reach the cold side.
        Depending on where we set the cold side, we check if phonon crossed that line"""

        small_offset = 10e-9
        if par.COLD_SIZE_POSITION == 'top':
            return par.LENGTH > self.y > 0
        if par.COLD_SIZE_POSITION == 'right':
            return (self.y < par.LENGTH - 1.1e-6) or (self.y > par.LENGTH - 1.1e-6 and self.x < par.WIDTH / 2.0 - small_offset)
        if par.COLD_SIZE_POSITION == 'top and right':
            return (self.y < 1.0e-6) or (1.0e-6 < self.y < par.LENGTH and self.x < par.WIDTH / 2.0 - small_offset)
        if par.COLD_SIZE_POSITION == 'top and bottom':
            return par.LENGTH > self.y > 0
        raise ValueError('Specified "cold_side" is not valid')

    def assign_initial_coordinates(self):
        """Assign initial coordinates at the hot side"""

        self.x = 0.49 * par.WIDTH * (2 * random() - 1)
        self.y = 1e-12
        self.z = 0.49 * par.THICKNESS * (2 * random() - 1)

        # In some lattices, phonons start from specific places:
        if par.HOLE_LATTICE_TYPE == 'serpentine':
            self.x = - par.WIDTH / 2 + 155e-9 / 2 + 0.4 * 155e-9 * (2 * random() - 1)
        if par.HOLE_LATTICE_TYPE == 'diode_with_wires':
            self.x = 0.4 * (par.PERIOD_X - par.RECTANGULAR_HOLE_SIDE_X) * (2 * random() - 1)
        if par.HOLE_LATTICE_TYPE == 'turn':
            self.x = 0.4 * (par.PERIOD_X * 5) * (2 * random() - 1)
        if par.HOLE_LATTICE_TYPE == 'turn90':
            self.x = 0.4 * par.WIDTH / 2.0 * (2 * random() - 1) - par.PERIOD_X * 3.5
        if par.HOLE_LATTICE_TYPE == 'directional_source':
            self.x = (2*random()-1)*50e-9
        if par.HOLE_LATTICE_TYPE == 'tesla valve forward':
            self.x = (2*random()-1)*10e-9 - 100e-9
        if par.HOLE_LATTICE_TYPE == 'tesla valve backward':
            self.x = (2*random()-1)*10e-9 - 100e-9
        if par.HOLE_LATTICE_TYPE == 'diamond_particle':
            self.x = 0.40 * par.THICKNESS * (2 * random() - 1)
        if par.HOLE_LATTICE_TYPE == 'cavity':
            self.y = 1550e-9


    def assign_angles(self):
        """Depending on angle distribution, assign angles"""
        if par.HOT_SIDE_ANGLE_DISTRIBUTION == 'random':
            self.theta = -pi/2 + pi*random()
            self.phi = asin(2*random() - 1)
        if par.HOT_SIDE_ANGLE_DISTRIBUTION == 'directional':
            self.theta = 0
            self.phi = -pi/2 + pi*random()
        if par.HOT_SIDE_ANGLE_DISTRIBUTION == 'lambert':
            self.theta = asin(2*random() - 1)
            self.phi = asin((asin(2*random() - 1))/(pi/2))
        if par.HOT_SIDE_ANGLE_DISTRIBUTION == 'uniform':
            self.theta = -pi + pi*2*random()
            self.phi = asin(2*random() - 1)


    def assign_frequency(self, material):
        """Assigning frequency with parobability according to Plankian distribution"""

        f_max = material.default_speed/(2*pi*hbar*material.default_speed/(2.82*k*par.T))# Frequency of the peak of the Plank distribution
        dos_max = 3*((2*pi*f_max)**2)/(2*(pi**2)*(material.default_speed**3))           # DOS for f_max in Debye apparoximation
        bose_einstein_max = 1/(exp((hbar*2*pi*f_max)/(k*par.T)) - 1)                    # Bose-Einstein distribution for f_max
        plank_distribution_max = dos_max*hbar*2*pi*f_max*bose_einstein_max              # Peak of the distribution (needed for normalization)

        while True:                                                                     # Until we obtain the frequency
            self.f = f_max*5*random()                                                   # Let's draw a random frequency in the 0 - 5*f_max range
            dos = 3*((2*pi*self.f)**2)/(2*(pi**2)*(material.default_speed**3))          # Calculate the DOS in Debye apparoximation
            bose_einstein = 1/(exp((hbar*2*pi*self.f)/(k*par.T))-1)                     # And the Bose-Einstein distribution
            plank_distribution = dos*hbar*2*pi*self.f*bose_einstein                     # Plank distribution

            # We take the normalized distribution at this frequency and draw a random number,
            # If the random number is lower than the distribution, we accept this frequency
            if random() < plank_distribution/plank_distribution_max and self.f < max(material.dispersion[:, 1]):
                break


    def assign_speed(self, material):
        """Calculate group velocity dw/dk according to the frequency and branch"""

        if self.polarization == 'TA' and self.f < max(material.dispersion[:,2]):        # If TA polarization and frequency is good
            j = abs((np.abs(material.dispersion[:, 2] - self.f)).argmin() - 1)           # Calculate index closest to f
            d_w = 2*pi*abs(material.dispersion[j+1, 2] - material.dispersion[j, 2])
        else:                                                                           # Otherwise it is LA polarization
            j = abs((np.abs(material.dispersion[:, 1] - self.f)).argmin() - 1)
            d_w = 2*pi*abs(material.dispersion[j+1, 1] - material.dispersion[j, 1])
        d_k = abs(material.dispersion[j+1,0] - material.dispersion[j,0])
        self.speed = d_w/d_k


    def assign_internal_scattering_time(self, material):
        """Determine relaxation time after which this phonon will undergo internal scattering"""

        if par.USE_GRAY_APPROXIMATION_MFP:
            self.time_of_internal_scattering = par.GRAY_APPROXIMATION_MFP / self.speed
        else:
            omega = 2*pi*self.f

            # Depending on the material we assign different relaxation times:
            if material.name == "Si":
                deb_temp = 152.0
                tau_impurity = 1/(2.95e-45 * (omega ** 4))
                tau_umklapp = 1/(0.95e-19 * (omega ** 2) * par.T * exp(-deb_temp / par.T))
                tau_internal = 1/((1/tau_impurity) + (1/tau_umklapp))

            elif material.name == "SiC":    # Ref. Joshi et al, JAP 88, 265 (2000)
                deb_temp = 1200
                tau_impurity = 1/(8.46e-45 * (omega ** 4.0))
                tau_umklapp = 1/(6.16e-20 * (omega ** 2.0) * par.T * exp(-deb_temp / par.T))
                tau_4p = 1/(6.9e-23 * (par.T ** 2) * (omega ** 2))
                tau_internal = 1/((1/tau_impurity) + (1/tau_umklapp) + (1/tau_4p))

            else:
                raise ValueError('Specified material does not exist in the database.')

        # Final relaxation time is determined with some randomization [PRB 94, 174303 (2016)]:
        self.time_of_internal_scattering = -log(random()) * tau_internal

    def move(self):
        """Move a phonon in one timestep and return new coordinates"""
        self.x, self.y, self.z = move.move(self)

    def correct_angle(self):
        """Check if angles are out of the [-pi:pi] range and return them back to this range"""
        if abs(self.theta) > pi:
            self.theta -= sign(self.theta)*2*pi
