"""This module provides phonon class which generates and moves a phonon"""

from math import pi, asin, exp, log
from random import random, choice
from numpy import sign
from scipy.constants import k, hbar
import numpy as np
import enum

from freepaths.config import cf
from freepaths.options import Distributions, Materials, Polarizations
import freepaths.move


class Phonon:
    """A phonon particle with various physical properties"""

    def __init__(self, material, polarization=None, phonon_number=None):
        """Initialize a phonon by assigning coordinates and other properties"""
        self.polarization = polarization
        self.phonon_number = phonon_number
        self.x = None
        self.y = None
        self.z = None
        self.f = None
        self.phi = None
        self.theta = None
        self.speed = None

        if polarization is None:
            self.assign_polarization()
        self.assign_coordinates()
        self.assign_angles()
        if phonon_number is None:
            self.assign_frequency(material)
        else:
            branch = polarization.value
            self.f = abs((material.dispersion[phonon_number+1, branch] + material.dispersion[phonon_number, branch]) / 2)
        self.assign_speed(material)
        self.assign_internal_scattering_time(material)

    @property
    def wavelength(self):
        """Calculate wavelength of the phonon"""
        return self.speed/self.f

    @property
    def is_in_system(self):
        """Checks if the phonon at this timestep did not reach the cold side.
        Depending on where user set cold sides, we check if phonon crossed that line"""
        is_inside_top = self.y < cf.length
        is_inside_bottom = self.y > 0
        is_inside_right = self.x < cf.width / 2.0
        is_inside_left = self.x > - cf.width / 2.0
        return ((not cf.cold_side_position_top or is_inside_top) and
                (not cf.cold_side_position_bottom or is_inside_bottom) and
                (not cf.cold_side_position_right or is_inside_right) and
                (not cf.cold_side_position_left or is_inside_left))

    def assign_polarization(self):
        """Assign branch of phonon dispersion"""
        self.polarization = choice([Polarizations.TA, Polarizations.TA, Polarizations.LA])

    def assign_coordinates(self):
        """Assign initial coordinates at the hot side"""
        self.x = cf.phonon_source_x + 0.49 * cf.phonon_source_width_x * (2 * random() - 1)
        self.y = cf.phonon_source_y + 0.49 * cf.phonon_source_width_y * (2 * random() - 1)
        self.z = 0.49 * cf.thickness * (2 * random() - 1)

    def assign_angles(self):
        """Depending on angle distribution, assign angles"""
        if cf.phonon_source_angle_distribution == Distributions.RANDOM_UP:
            self.theta = -pi/2 + pi*random()
            self.phi = asin(2*random() - 1)
        if cf.phonon_source_angle_distribution == Distributions.RANDOM_DOWN:
            rand_sign = sign((2*random() - 1))
            self.theta = rand_sign*(pi/2 + pi/2*random())
            self.phi = asin(2*random() - 1)
        if cf.phonon_source_angle_distribution == Distributions.RANDOM_RIGHT:
            self.theta = pi*random()
            self.phi = asin(2*random() - 1)
        if cf.phonon_source_angle_distribution == Distributions.RANDOM_LEFT:
            self.theta = - pi*random()
            self.phi = asin(2*random() - 1)
        if cf.phonon_source_angle_distribution == Distributions.DIRECTIONAL:
            self.theta = 0
            self.phi = -pi/2 + pi*random()
        if cf.phonon_source_angle_distribution == Distributions.LAMBERT:
            self.theta = asin(2*random() - 1)
            self.phi = asin((asin(2*random() - 1))/(pi/2))
        if cf.phonon_source_angle_distribution == Distributions.UNIFORM:
            self.theta = -pi + 2*pi*random()
            self.phi = asin(2*random() - 1)

    def assign_frequency(self, material):
        """Assigning frequency with probability according to Planckian distribution"""

        # Frequency of the peak of the Plank distribution:
        f_max = material.default_speed/(2*pi*hbar*material.default_speed/(2.82*k*cf.temp))

        # Density of states for f_max in Debye apparoximation:
        dos_max = 3*((2*pi*f_max)**2)/(2*(pi**2)*(material.default_speed**3))

        # Bose-Einstein distribution for f_max:
        bose_einstein_max = 1/(exp((hbar*2*pi*f_max)/(k*cf.temp)) - 1)

        # Peak of the distribution (needed for normalization)
        plank_distribution_max = dos_max*hbar*2*pi*f_max*bose_einstein_max

        while True:                                                             # Until we obtain the frequency
            self.f = f_max*5*random()                                           # Let's draw a random frequency in the 0 - 5*f_max range
            dos = 3*((2*pi*self.f)**2)/(2*(pi**2)*(material.default_speed**3))  # Calculate the DOS in Debye apparoximation
            bose_einstein = 1/(exp((hbar*2*pi*self.f)/(k*cf.temp))-1)           # And the Bose-Einstein distribution
            plank_distribution = dos*hbar*2*pi*self.f*bose_einstein             # Plank distribution

            # We take the normalized distribution at this frequency and draw a random number,
            # If the random number is lower than the distribution, we accept this frequency
            if random() < plank_distribution/plank_distribution_max and self.f < max(material.dispersion[:, 1]):
                break

    def assign_speed(self, material):
        """Calculate group velocity dw/dk according to the frequency and polarization"""
        if self.polarization == Polarizations.TA and self.f < max(material.dispersion[:,2]):
            point_num = abs((np.abs(material.dispersion[:, 2] - self.f)).argmin() - 1)
            d_w = 2*pi*abs(material.dispersion[point_num+1, 2] - material.dispersion[point_num, 2])
        else:
            point_num = abs((np.abs(material.dispersion[:, 1] - self.f)).argmin() - 1)
            d_w = 2*pi*abs(material.dispersion[point_num+1, 1] - material.dispersion[point_num, 1])
        d_k = abs(material.dispersion[point_num+1, 0] - material.dispersion[point_num, 0])
        self.speed = d_w/d_k

    def assign_internal_scattering_time(self, material):
        """Determine relaxation time after which this phonon will undergo internal scattering"""

        if cf.use_gray_approximation_mfp:
            self.time_of_internal_scattering = cf.gray_approximation_mfp / self.speed
        else:
            omega = 2*pi*self.f

            # Depending on the material we assign different relaxation times:
            if material.name == Materials.Si:
                deb_temp = 152.0
                tau_impurity = 1/(2.95e-45 * (omega ** 4))
                tau_umklapp = 1/(0.95e-19 * (omega ** 2) * cf.temp * exp(-deb_temp / cf.temp))
                tau_internal = 1/((1/tau_impurity) + (1/tau_umklapp))

            elif material.name == Materials.SiC:    # Ref. Joshi et al, JAP 88, 265 (2000)
                deb_temp = 1200
                tau_impurity = 1/(8.46e-45 * (omega ** 4))
                tau_umklapp = 1/(6.16e-20 * (omega ** 2) * cf.temp * exp(-deb_temp / cf.temp))
                tau_4p = 1/(6.9e-23 * (cf.temp ** 2) * (omega ** 2))
                tau_internal = 1/((1/tau_impurity) + (1/tau_umklapp) + (1/tau_4p))

            else:
                raise ValueError('Specified material does not exist in the database.')

            # Final relaxation time is determined with some randomization [PRB 94, 174303 (2016)]:
            self.time_of_internal_scattering = -log(random()) * tau_internal

    def move(self):
        """Move a phonon in one timestep and return new coordinates"""
        self.x, self.y, self.z = freepaths.move.move(self, cf.timestep)

    def correct_angle(self):
        """Check if angles are out of the [-pi:pi] range and return them back to this range"""
        if abs(self.theta) > pi:
            self.theta -= sign(self.theta)*2*pi
