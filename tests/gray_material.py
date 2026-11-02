"""
A test-only "gray" material: three identical linear (Debye) branches with one group velocity,
and one constant relaxation time, all inelastic. Every mode then has the same mean free path
Lambda = v*tau, and the bulk thermal conductivity is exactly kappa = C v^2 tau / 3, with C the
heat capacity of the tabulated branches.

install() registers it under the name "Gray" in freepaths.materials. An input file calls it,
so that it runs in every worker process too (each process executes the input file).
"""

import numpy as np

import freepaths.materials as materials


VELOCITY = 6000.0          # [m/s]
RELAXATION_TIME = 10e-12   # [s], so Lambda = 60 nm
K_MAX = 1.2e10             # [1/m]


class Gray(materials.Si):
    """Si with linear branches and a constant, purely inelastic relaxation time"""

    def assign_phonon_dispersion(self, num_points):
        k = np.linspace(0, K_MAX, num_points)
        f = VELOCITY * k / (2 * np.pi)
        self.dispersion = np.column_stack((k, f, f, f))

    def phonon_scattering_rates(self, omega):
        return 1 / RELAXATION_TIME, 0.0

    def phonon_relaxation_time(self, omega):
        return RELAXATION_TIME


def install():
    original = materials.get_media_class
    if getattr(original, "with_gray", False):
        return

    def get_media_class(name):
        return Gray if name == "Gray" else original(name)
    get_media_class.with_gray = True
    materials.get_media_class = get_media_class


def bulk_conductivity(temp):
    """kappa = C v^2 tau / 3 with C from the same tabulated branches the simulation samples"""
    return Gray(temp).dispersion_heat_capacity * VELOCITY**2 * RELAXATION_TIME / 3
