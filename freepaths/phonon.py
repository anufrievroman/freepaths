"""This module provides phonon class which generates and moves a phonon"""

from math import pi, asin, exp, log
from random import random, choice, randint
from numpy import sign
from scipy.constants import k, hbar
import numpy as np
import enum

from freepaths.config import cf
import freepaths.move


class Phonon:
    """A phonon particle with various physical properties"""

    def __init__(self, material, branch_number=None, phonon_number=None):
        """Initialize a phonon by assigning initial properties"""
        self.branch_number = branch_number
        self.phonon_number = phonon_number
        self.x = None
        self.y = None
        self.z = None
        self.f = None
        self.phi = None
        self.theta = None
        self.speed = None
        self.first_timestep = randint(0, cf.number_of_virtual_timesteps)

        # Assign initial properties of the phonon:
        if self.branch_number is None:
            self.branch_number = choice(range(3))
        self.f_max = max(material.dispersion[:, self.branch_number + 1])

        # Assign initial coordinates but ensure that it's not inside a hole:
        source = choice(cf.phonon_sources)
        while True:
            self.x, self.y, self.z = source.generate_coordinates()
            is_in_hole = any(hole.is_inside(self.x, self.y, None, cf) for hole in cf.holes)
            if not is_in_hole:
                break

        # Assign initial angles:
        self.theta, self.phi = source.generate_angles()
        self.correct_angle()
        if cf.is_two_dimensional_material:
            self.phi = 0.0
            self.z = 0.0

        # Frequency is assigned based on Plankian distribution in case MFP sampling mode:
        if phonon_number is None:
            self.assign_frequency(material)

        # Otherwise, frequency is just assigned depending on the phonon number:
        else:
            f_upper = abs(material.dispersion[phonon_number + 1, self.branch_number + 1])
            f_lower = abs(material.dispersion[phonon_number, self.branch_number +1])
            self.f = (f_upper + f_lower) / 2

        self.assign_speed(material)
        self.assign_internal_scattering_time(material)

    @property
    def wavelength(self):
        """Calculate wavelength of the phonon"""
        return self.speed / self.f

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

    def assign_frequency(self, material):
        """Assigning frequency with probability according to Planckian distribution"""

        # Frequency of the peak of the Plank distribution:
        f_peak = material.default_speed/(2*pi*hbar*material.default_speed/(2.82*k*cf.temp))

        # Density of states for f_max in Debye approximation:
        dos_max = 3*((2*pi*f_peak)**2)/(2*(pi**2)*(material.default_speed**3))

        # Bose-Einstein distribution for f_max:
        bose_einstein_max = 1/(exp((hbar*2*pi*f_peak)/(k*cf.temp)) - 1)

        # Peak of the distribution (needed for normalization)
        plank_distribution_max = dos_max*hbar*2*pi*f_peak*bose_einstein_max

        while True:                                                 # Until we obtain the frequency
            self.f = f_peak*5*random()                              # Let's draw a random frequency in the 0 - 5*f_peak range
            dos = 3*((2*pi*self.f)**2)/(2*(pi**2)*(material.default_speed**3))  # Calculate the DOS in Debye approximation
            bose_einstein = 1/(exp((hbar*2*pi*self.f)/(k*cf.temp))-1)           # And the Bose-Einstein distribution
            plank_distribution = dos*hbar*2*pi*self.f*bose_einstein             # Plank distribution

            # We take the normalized distribution at this frequency and draw a random number,
            # If the random number is lower than the distribution, we accept this frequency
            if random() < plank_distribution/plank_distribution_max and self.f < self.f_max:
                break

    def assign_speed(self, material):
        """Calculate group velocity dw/dk according to the frequency and polarization"""
        point_num = abs((np.abs(material.dispersion[:, self.branch_number + 1] - self.f)).argmin() - 1)
        d_w = 2*pi*abs(material.dispersion[point_num + 1, self.branch_number + 1]
                       - material.dispersion[point_num, self.branch_number + 1])
        d_k = abs(material.dispersion[point_num + 1, 0] - material.dispersion[point_num, 0])
        self.speed = d_w / d_k

    def assign_internal_scattering_time(self, material):
        """Determine relaxation time after which this phonon will undergo internal scattering"""
        if cf.use_gray_approximation_mfp:
            self.time_of_internal_scattering = cf.gray_approximation_mfp / self.speed
        else:
            # Relaxation time is assigned with some randomization [PRB 94, 174303 (2016)]:
            omega = 2 * pi * self.f
            self.time_of_internal_scattering = -log(random()) * material.relaxation_time(omega)

    def move(self):
        """Move a phonon in one timestep and return new coordinates"""
        self.x, self.y, self.z = freepaths.move.move(self, cf.timestep)

    def correct_angle(self):
        """Check if angles are out of the [-pi:pi] range and return them back to this range"""
        if abs(self.theta) > pi:
            self.theta -= sign(self.theta)*2*pi
