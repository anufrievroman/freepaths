"""This module provides phonon class which generates and moves a phonon"""

from math import pi, exp, log
from random import random, choice
from numpy import sign
from scipy.constants import k, hbar, h
import numpy as np

from freepaths.config import cf
from freepaths.options import ParticleType, SimulationMode
from freepaths.particle import Particle
import freepaths.move


class Phonon(Particle):
    """A phonon particle with various physical properties"""

    def __init__(self, material, mode, branch_number=None, phonon_number=None):
        """Initialize a phonon by assigning initial properties"""
        super().__init__()
        # Assign particle type
        self.type = ParticleType.PHONON
        self.mode = mode

        self.branch_number = branch_number
        self.phonon_number = phonon_number

        # Assign initial coordinates but ensure that it's not inside a hole:
        source = choice(cf.particles_sources)
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

        # Assign branch and frequency. Each branch below sets self.branch_number itself,
        # except MFP sampling, where the caller (main_mfp_sampling.py) already fixed it:
        if mode is SimulationMode.PHONON_MFP_SAMPLING:
            self.assign_mfp_sampling_frequency(material)
        else:
            self.assign_frequency(material)

        # draw_from_table already sets the speed from the drawn dispersion bin:
        if self.speed is None:
            self.assign_speed(material)
        self.assign_internal_scattering_time(material)

    @property
    def wavelength(self):
        """Calculate wavelength of the phonon"""
        return self.speed / self.f

    @property
    def energy(self):
        """Calculate energy of the phonon"""
        return self.f * h

    def assign_mfp_sampling_frequency(self, material):
        """Set frequency deterministically from the phonon index sweeping the dispersion grid (MFP sampling mode; branch_number is fixed by the caller)."""
        f_upper = abs(material.dispersion[self.phonon_number + 1, self.branch_number + 1])
        f_lower = abs(material.dispersion[self.phonon_number, self.branch_number + 1])
        self.f = (f_upper + f_lower) / 2

    def assign_frequency(self, material):
        """
        Assign branch and frequency for phonon tracing mode. Uses the pre-tabulated,
        dispersion-weighted cumulative probabilities by default (see Material.assign_phonon_sampling_tables); falls back to a
        branch-blind Debye-Planck rejection sampling under the gray approximation or
        when SAMPLE_FROM_DISPERSION is off, both of which treat the 3 branches as
        identical.
        """
        if cf.sample_from_dispersion and not cf.use_gray_approximation_mfp:
            # Particles are equal-energy bundles of deviational energy, so the source
            # emits with the heat-capacity weight C(w), consistent with rethermalize():
            self.draw_from_table(material, material.emission_branch_probabilities, material.emission_frequency_probabilities)
        else:
            self.assign_debye_frequency(material)

    def assign_debye_frequency(self, material):
        """
        Legacy/gray path: pick a branch uniformly at random (the lumped Debye model
        below treats all 3 branches as identical) and draw a frequency via rejection
        sampling from the Planck distribution, weighted by the Debye density of states.
        """
        # Branch is not weighted here since the Debye model below is branch-blind
        self.branch_number = choice(range(3))
        self.f_max = max(material.dispersion[:, self.branch_number + 1])

        # Frequency of the peak of the Plank distribution:
        f_peak = material.default_speed/(2*pi*hbar*material.default_speed/(2.82*k*cf.temp))

        # Density of states for f_max in Debye approximation:
        dos_max = 3*((2*pi*f_peak)**2)/(2*(pi**2)*(material.default_speed**3))

        # Bose-Einstein distribution for f_max:
        bose_einstein_max = 1/(exp((hbar*2*pi*f_peak)/(k*cf.temp)) - 1)

        # Peak of the distribution (needed for normalization)
        plank_distribution_max = dos_max*hbar*2*pi*f_peak*bose_einstein_max

        while True:
            # Draw a candidate frequency in the 0 - 5*f_peak range
            self.f = f_peak*5*random()
            # Density of states in Debye approximation
            dos = 3*((2*pi*self.f)**2)/(2*(pi**2)*(material.default_speed**3))
            # Bose-Einstein distribution at this frequency
            bose_einstein = 1/(exp((hbar*2*pi*self.f)/(k*cf.temp))-1)
            # Planck distribution at this frequency
            plank_distribution = dos*hbar*2*pi*self.f*bose_einstein

            # Accept the candidate with probability proportional to the normalized distribution
            if random() < plank_distribution/plank_distribution_max and self.f < self.f_max:
                break

    def draw_from_table(self, material, branch_probabilities, frequency_probabilities):
        """
        Draw a branch and frequency by inverse sampling of one of the material's
        pre-tabulated cumulative probabilities (see Material.assign_phonon_sampling_tables).
        """
        # Pick the branch via inverse sampling of the cumulative branch probabilities
        self.branch_number = int(np.searchsorted(branch_probabilities, random()))
        # Upper frequency bound for this branch, used elsewhere to reject invalid draws
        self.f_max = max(material.dispersion[:, self.branch_number + 1])
        freqs = material.frequencies_table[self.branch_number]
        # Pick the frequency bin via inverse sampling within the chosen branch
        index = min(np.searchsorted(frequency_probabilities[self.branch_number], random()), len(freqs) - 1)
        self.f = freqs[index]
        # Group velocity of exactly the drawn dispersion bin:
        self.speed = material.group_velocity_table[self.branch_number][index]

    def rethermalize(self, material):
        """
        Re-draw branch and frequency at an internal scattering event from the
        collision-rate-weighted distribution C(f)/tau(f) (see Material.assign_phonon_sampling_tables), following
        Peraud & Hadjiconstantinou, PRB 84, 205331 (2011). The 1/tau weight is
        essential: phonons drawn from this distribution spend a time tau(f) at
        each frequency, so the steady-state population spectrum remains equal
        to the emitted spectrum and does not drift.
        """
        # Draw the new branch and frequency from the scattering tables;
        # the group velocity of the drawn bin is set by draw_from_table
        self.draw_from_table(material, material.scattering_branch_probabilities, material.scattering_frequency_probabilities)

    def assign_speed(self, material):
        """Set group velocity dw/dk at the phonon's frequency and branch."""
        self.speed = material.group_velocity(self.branch_number, self.f)

    def assign_internal_scattering_time(self, material):
        """Determine relaxation time after which this phonon will undergo internal scattering"""
        if cf.use_gray_approximation_mfp:
            self.time_of_internal_scattering = cf.gray_approximation_mfp / self.speed
        else:
            # Relaxation time is assigned with some randomization [PRB 94, 174303 (2016)]:
            omega = 2 * pi * self.f
            self.time_of_internal_scattering = -log(random()) * material.phonon_relaxation_time(omega)

    def move(self):
        """Move a phonon in one timestep and return new coordinates"""
        self.x, self.y, self.z = freepaths.move.move(self, cf.timestep)

    def correct_angle(self):
        """Check if angles are out of the [-pi:pi] range and return them back to this range"""
        if abs(self.theta) > pi:
            self.theta -= sign(self.theta)*2*pi
