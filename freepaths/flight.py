"""Module that provides particle flight characteristics"""

from math import cos
import numpy as np

from freepaths.config import cf


class Path:
    """Particle path coordinates"""

    def __init__(self, x_init, y_init, z_init):
        """Initialize a path"""
        self.x = [x_init]
        self.y = [y_init]
        self.z = [z_init]

    def add_point(self, x_new, y_new, z_new):
        """Add a point to the particle path"""
        self.x.append(x_new)
        self.y.append(y_new)
        self.z.append(z_new)

    @property
    def number_of_path_points(self):
        """Calculate number of points in the path"""
        return len(self.x)


class Flight:
    """Particle flight characteristics"""

    def __init__(self, particle):
        """Initialize a particle
        flight"""
        self.particle = particle
        self.initial_frequency = self.particle.f
        self.initial_theta = self.particle.theta
        self.path = Path(self.particle.x, self.particle.y, self.particle.z)
        self.exit_theta = 0.0
        self.exit_frequency = 0.0
        self.free_path = 0.0
        self.travel_time = 0.0
        self.time_since_previous_scattering = 0.0
        self.free_paths = []
        self._mfp_sum = 0.0
        self._mfp_count = 0
        self.hole_diff_scattering_angles = []
        self.hole_spec_scattering_angles = []
        self.thermal_conductivity = 0.0
    @property
    def mean_free_path(self):
        """Mean value of all free flights"""
        return self._mfp_sum / self._mfp_count if self._mfp_count else 0

    @property
    def mean_free_path_sem(self):
        """Standard error of the mean free path (MFP sampling mode only). Assumes
        individual free-path segments are exponentially distributed (the natural
        assumption for Poisson-process intrinsic scattering: coefficient of
        variation = 1, i.e. std = mean), so SEM = mean_free_path / sqrt(n).

        An earlier version instead measured the spread empirically from this one
        phonon's own segments (a running sum of squares, no per-segment storage
        needed) - but that returns exactly 0 for any mode with fewer than 2
        segments, silently zeroing out the modes that reach the cold side in a
        single, boundary-limited hop - precisely the near-ballistic population
        with the *most* run-to-run variability. That version underestimated the
        true spread by ~5x (rerunning one config - Nanowire_Si_MFP_L5.5um, N=200,
        cap=1000 - 10 times gave empirical kappa std = 1.08 W/mK vs. its ~0.21
        W/mK average prediction). This exponential-CV version was validated the
        same way (9 reruns): empirical std 0.91 vs. 0.94 W/mK average prediction
        (session notes, CLAUDE.md, 2026-07-20)."""
        if self._mfp_count < 1:
            return 0.0
        return self.mean_free_path / (self._mfp_count ** 0.5)

    def add_point_to_path(self):
        """Add a scattering point to the path"""
        self.path.add_point(self.particle.x, self.particle.y, self.particle.z)

    def save_free_paths(self):
        """Save current free path to the list of free paths"""
        self._mfp_sum += self.free_path
        self._mfp_count += 1
        if not cf.low_memory_usage:
            self.free_paths.append(self.free_path)

    def save_hole_diff_scattering_angle(self, angle):
        """Save angle of diffuse scattering from the hole"""
        self.hole_diff_scattering_angles.append(angle)

    def save_hole_spec_scattering_angle(self, angle):
        """Save angle of specular scattering from the hole"""
        self.hole_spec_scattering_angles.append(angle)

    def restart(self):
        """Restart the flight after a scattering event"""
        self.time_since_previous_scattering = 0.0
        self.free_path = 0.0

    def finish(self, step, timestep, reached_cold_side=False):
        """Finish the flight and record final state. travel_time is only set when
        the particle genuinely reached the cold side - it means "time to cold side",
        so leaving it at 0 for other exits (e.g. the MFP-sampling scattering-event
        cap) keeps that meaning intact instead of recording a meaningless duration."""
        self.exit_theta = self.particle.theta
        self.exit_frequency = self.particle.f
        if reached_cold_side:
            self.travel_time = step * timestep

    def add_step(self, timestep):
        """Increase parameters of the flight by length of one step"""
        step_length = self.particle.speed * timestep
        self.free_path += step_length
        self.time_since_previous_scattering += timestep

    def reset_travel_time(self):
        """Reset the travel time of a particle"""
        self.travel_time = 0.0
