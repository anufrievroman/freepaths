"""This module provides phonon class which generates and moves a phonon"""

from math import pi, exp, log, sin, cos, sqrt, asin, atan2, expm1
from random import random, choice
from numpy import sign
from scipy.constants import h, hbar, k as k_B
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
        self.grain_size = self._draw_grain_size()
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
        Assign branch and frequency for phonon tracing mode by inverse sampling of the
        pre-tabulated, dispersion-weighted cumulative probabilities
        (see Material.assign_phonon_sampling_tables). Particles are equal-energy bundles
        of deviational energy, so the source emits with the heat-capacity weight C(w),
        consistent with rethermalize().
        """
        self.draw_from_table(material, material.emission_branch_probabilities, material.emission_frequency_probabilities)

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

    def _draw_grain_size(self):
        """Draw a grain diameter from lognormal(GRAIN_SIZE, GRAIN_SIZE_STD), or return None."""
        if cf.grain_size is None:
            return None
        if cf.grain_size_std > 0:
            s2 = np.log(1 + (cf.grain_size_std / cf.grain_size) ** 2)
            mu = np.log(cf.grain_size) - s2 / 2
            return np.random.lognormal(mu, np.sqrt(s2))
        return cf.grain_size

    def _grain_boundary_rate(self):
        """Grain boundary scattering rate with Soffer-type low-frequency transparency.
        Rate = (v_g / d) * (1 - exp(-4 sigma_GB^2 omega^2 / v_g^2)).
        Vanishes at low omega (long wavelength passes through); approaches v_g/d at high omega."""
        if self.grain_size is None:
            return 0.0
        omega = 2 * pi * self.f
        soffer = exp(-4 * cf.grain_roughness ** 2 * omega ** 2 / self.speed ** 2)
        return (self.speed / self.grain_size) * (1 - soffer)

    def classify_internal_event(self, material):
        """
        Attribute the current internal scattering event to one of three channels.
        The event time is drawn from the total rate (a sum of independent Poisson
        processes), so the event lands in each channel with probability equal to its
        share of the total rate:
        - "normal": momentum-conserving anharmonic (Normal) three-phonon process. Only
          active in hydrodynamic mode; redraws the mode and a drift-biased direction
          (drift_scatter), conserving the local crystal momentum.
        - "inelastic": momentum-destroying anharmonic (Umklapp / 4-phonon) process;
          rethermalizes the mode (isotropic redraw).
        - "elastic": point-defect / impurity / grain-boundary process; conserves the
          mode and only randomizes the direction (already done by internal_scattering).
        Returns one of the strings above. When hydrodynamics is off the normal rate is
        zero, so this reduces exactly to the previous inelastic/elastic split.
        """
        inelastic_rate, elastic_rate = material.phonon_scattering_rates(2 * pi * self.f)
        elastic_rate += self._grain_boundary_rate()
        normal_rate = material.phonon_normal_rate(2 * pi * self.f) if cf.phonon_hydrodynamic else 0.0
        total_rate = normal_rate + inelastic_rate + elastic_rate
        if total_rate == 0.0:
            return "inelastic"
        threshold = random() * total_rate
        if threshold < normal_rate:
            return "normal"
        if threshold < normal_rate + inelastic_rate:
            return "inelastic"
        return "elastic"

    def drift_scatter(self, material, drift_velocity):
        """
        Momentum-conserving Normal scattering event. The equal-energy bundle keeps its
        deviational energy, so only the mode label and direction change: the new mode is
        drawn from the collision-rate-weighted table (as in rethermalize, since the
        Normal rate shares the C/tau frequency dependence of the anharmonic rate), and
        the new direction is drawn from the local drifting distribution so that the
        outgoing population carries the local crystal momentum rather than relaxing it.
        """
        self.draw_from_table(material, material.scattering_branch_probabilities, material.scattering_frequency_probabilities)
        # Diagnostic control: strip the drift bias so N redraws isotropically (momentum-
        # destroying), while everything else (rate, mode redraw, momentum recording) stays:
        if cf.hydrodynamic_normal_resistive:
            drift_velocity = (0.0, 0.0, 0.0)
        self._draw_drift_direction(material, drift_velocity)

    def _draw_drift_direction(self, material, drift_velocity):
        """
        Draw a new propagation direction from the linearized displaced-Bose-Einstein
        angular distribution P(d_hat) ~ 1 + a * (d_hat . u_hat), where the bias
        a = (1 + n0) * hbar * k * |u| / (k_B T) aligns the outgoing momentum with the
        local drift u (clamped to [0, 1] to keep P >= 0 in the small-drift regime). With
        |u| = 0 (bootstrap pass) this reduces to an isotropic redraw, i.e. momentum-
        neutral. Sampled by rejection; the polar angle is measured from u_hat.
        """
        u_x, u_y, u_z = drift_velocity
        u_mag = sqrt(u_x**2 + u_y**2 + u_z**2)

        # Bootstrap / no drift: isotropic redraw (momentum-neutral).
        if u_mag == 0.0:
            self.theta = -pi + random() * 2 * pi
            self.phi = 0.0 if cf.is_two_dimensional_material else asin(2 * random() - 1)
            return

        omega = 2 * pi * self.f
        k = material.wavevector(self.branch_number, self.f)
        n0 = 1.0 / expm1(hbar * omega / (k_B * material.temp))
        a = (1.0 + n0) * hbar * k * u_mag / (k_B * material.temp)
        a = min(a, 1.0)

        if cf.is_two_dimensional_material:
            # In-plane: sample the angle delta from u_hat with P(delta) ~ 1 + a cos(delta).
            theta_u = atan2(u_x, u_y)
            while True:
                delta = -pi + random() * 2 * pi
                if random() * (1.0 + a) <= 1.0 + a * cos(delta):
                    break
            self.theta = theta_u + delta
            self.phi = 0.0
            self.correct_angle()
            return

        # 3D: sample cos(psi) from P ~ 1 + a cos(psi), azimuth uniform about u_hat.
        while True:
            cos_psi = 2 * random() - 1
            if random() * (1.0 + a) <= 1.0 + a * cos_psi:
                break
        sin_psi = sqrt(max(0.0, 1.0 - cos_psi**2))
        azimuth = 2 * pi * random()

        # Build an orthonormal frame (e1, e2) perpendicular to u_hat:
        u_hat = np.array([u_x, u_y, u_z]) / u_mag
        helper = np.array([1.0, 0.0, 0.0]) if abs(u_hat[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        e1 = helper - helper.dot(u_hat) * u_hat
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(u_hat, e1)
        d = cos_psi * u_hat + sin_psi * (cos(azimuth) * e1 + sin(azimuth) * e2)

        # Convert direction d to (theta, phi) matching the move convention
        # d = (sin(theta) |cos(phi)|, cos(theta) |cos(phi)|, sin(phi)):
        self.phi = asin(max(-1.0, min(1.0, d[2])))
        horizontal = sqrt(max(0.0, 1.0 - d[2]**2))
        self.theta = atan2(d[0], d[1]) if horizontal > 1e-12 else 0.0
        self.correct_angle()

    def rethermalize(self, material):
        """
        Re-draw branch and frequency at an inelastic internal scattering event from
        the collision-rate-weighted distribution C(f)/tau_inelastic(f) (see
        Material.assign_phonon_sampling_tables), following Peraud &
        Hadjiconstantinou, PRB 84, 205331 (2011). The 1/tau_inelastic weight is
        essential: with events fired at the total rate but redrawn only when the
        channel draw lands on the inelastic one, a phonon waits on average exactly
        tau_inelastic(f) at each drawn frequency, so the steady-state population
        spectrum remains equal to the emitted spectrum and does not drift.
        """
        # Draw the new branch and frequency from the scattering tables;
        # the group velocity of the drawn bin is set by draw_from_table
        self.draw_from_table(material, material.scattering_branch_probabilities, material.scattering_frequency_probabilities)

    def assign_speed(self, material):
        """Set group velocity dw/dk at the phonon's frequency and branch."""
        self.speed = material.group_velocity(self.branch_number, self.f)

    def assign_internal_scattering_time(self, material):
        """Determine relaxation time after which this phonon will undergo internal scattering"""
        # Relaxation time is assigned with some randomization [PRB 94, 174303 (2016)]:
        omega = 2 * pi * self.f
        total_rate = 1 / material.phonon_relaxation_time(omega) + self._grain_boundary_rate()
        # In hydrodynamic mode the momentum-conserving Normal events also scatter the
        # phonon (they redirect it, conserving momentum), so they must fire as part of
        # the Poisson clock. Excluded from phonon_relaxation_time (non-resistive), so
        # they are added here explicitly:
        if cf.phonon_hydrodynamic:
            total_rate += material.phonon_normal_rate(omega)
        self.time_of_internal_scattering = -log(random()) / total_rate

    def move(self):
        """Move a phonon in one timestep and return new coordinates"""
        self.x, self.y, self.z = freepaths.move.move(self, cf.timestep)

    def correct_angle(self):
        """Check if angles are out of the [-pi:pi] range and return them back to this range"""
        if abs(self.theta) > pi:
            self.theta -= sign(self.theta)*2*pi
