"""
A test-only non-gray material without mode mixing: three identical linear branches with one
group velocity, and a purely ELASTIC relaxation time that takes one value below K_SPLIT and a
shorter one above it. A phonon never changes its mode inside the medium, so every mode diffuses
on its own and the bulk thermal conductivity is exactly kappa = sum_i C_i v^2 tau_i / 3.

Unlike the gray material, it tells apart the two ways a hot side can treat the phonons that
come back to it: absorbing them and emitting new ones (correct, the arithmetic average above),
or re-emitting them in the same mode (every mode then carries heat in proportion to its
emission rate only, which gives the harmonic-type average of harmonic_conductivity()).

install() registers it under the name "TwoTau" in freepaths.materials, like gray_material.
"""

import numpy as np
from scipy.constants import hbar, k as k_B

import freepaths.materials as materials


VELOCITY = 6000.0        # [m/s]
K_MAX = 1.2e10           # [1/m]
K_SPLIT = 0.6 * K_MAX    # [1/m]
TAU_LONG = 10e-12        # [s] below K_SPLIT, Lambda = 60 nm
TAU_SHORT = 2e-12        # [s] above K_SPLIT, Lambda = 12 nm


def relaxation_time(omega):
    return TAU_LONG if omega / VELOCITY < K_SPLIT else TAU_SHORT


class TwoTau(materials.Si):
    """Si with linear branches and a two-valued, elastic relaxation time"""

    def assign_phonon_dispersion(self, num_points):
        k = np.linspace(0, K_MAX, num_points)
        f = VELOCITY * k / (2 * np.pi)
        self.dispersion = np.column_stack((k, f, f, f))

    def phonon_scattering_rates(self, omega):
        # A negligible inelastic rate, only to keep the rethermalization table defined:
        return 1.0, 1 / relaxation_time(omega)

    def phonon_relaxation_time(self, omega):
        return relaxation_time(omega)


def install():
    original = materials.get_media_class
    if getattr(original, "with_two_tau", False):
        return

    def get_media_class(name):
        return TwoTau if name == "TwoTau" else original(name)
    get_media_class.with_two_tau = True
    materials.get_media_class = get_media_class


def _mode_heat_capacities_and_times(temp):
    k = TwoTau(temp).dispersion[:, 0]
    k_mid, d_k = (k[1:] + k[:-1]) / 2, np.diff(k)
    omega = VELOCITY * k_mid
    x = hbar * omega / (k_B * temp)
    heat_capacity = 3 * k_B * x**2 * np.exp(x) / np.expm1(x)**2 * k_mid**2 * d_k / (2 * np.pi**2)
    return heat_capacity, np.array([relaxation_time(w) for w in omega])


def bulk_conductivity(temp):
    """kappa = sum C_i v^2 tau_i / 3 over the tabulated modes"""
    heat_capacity, tau = _mode_heat_capacities_and_times(temp)
    return np.sum(heat_capacity * tau) * VELOCITY**2 / 3


def harmonic_conductivity(temp):
    """What a hot side that re-emits the returning phonons in the same mode would give"""
    heat_capacity, tau = _mode_heat_capacities_and_times(temp)
    return np.sum(heat_capacity)**2 / np.sum(heat_capacity / tau) * VELOCITY**2 / 3
