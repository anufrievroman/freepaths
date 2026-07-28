"""This module provides electron class which generates and moves an electron"""

from random import choice, random
from scipy.constants import h, k, hbar, electron_volt, elementary_charge, epsilon_0
from freepaths.config import cf
from freepaths.particle import Particle
from freepaths.options import ParticleType

import numpy as np


def _carrier_dos_mass(material):
    """DOS effective mass of the active carrier (electron or hole)."""
    return material.effective_electron_dos_mass if cf.is_carrier_electron else material.effective_hole_dos_mass


def depletion_width(material):
    """Surface-depletion (dead-layer) width for carriers, from the abrupt-junction
    depletion approximation W = sqrt(2*eps_s*phi_s/(e*N_D)), with phi_s = SURFACE_POTENTIAL
    the surface band bending (eV, used as volts), N_D = DOPING_CONCENTRATION the doping,
    and eps_s the static permittivity of the material. Carriers are excluded within W of
    every free surface (Fermi-level pinning by surface states / native oxide). Returns 0
    if there is no doping or no surface potential. Reference: any semiconductor-device text
    (e.g. Sze), surface-depletion / dead-layer model for nanostructure conductivity."""
    N_I = cf.doping_concentration
    phi = cf.surface_potential
    eps_r = getattr(material, 'dielectric_constant', None)
    if not N_I or phi <= 0 or not eps_r:
        return 0.0
    eps = eps_r * epsilon_0
    return np.sqrt(2.0 * eps * phi / (elementary_charge * N_I))


def electron_mfp(energy, material):
    """Energy- and doping-dependent carrier mean free path via Matthiessen's rule:
        1/Lambda_tot(E) = 1/Lambda_ac + 1/Lambda_ii(E, N_I).

    Lambda_ac = ELECTRON_MFP is the crystal (acoustic-phonon-limited) MFP, taken
    energy-independent (for acoustic deformation-potential scattering tau ~ E^-1/2
    cancels v ~ E^1/2, so Lambda = v*tau ~ E^0). Lambda_ii is the ionized-impurity
    contribution from the Brooks-Herring model: screened Coulomb scattering off the
    DOPING_CONCENTRATION ionized dopants, with tau_ii ~ E^{3/2} in the weak-screening
    limit. Ionized impurities thus scatter LOW-energy carriers most, shortening the
    MFP with doping and making scattering energy-selective, which raises the Seebeck
    coefficient. DOPING_CONCENTRATION = 0 recovers the crystal MFP (optionally with
    the legacy phenomenological exponent ELECTRON_MFP_ENERGY_EXPONENT).

    Applies identically to electrons and holes: the carrier's own DOS effective mass
    is used (chosen by IS_CARRIER_ELECTRON) and N_I is the density of ionized donors
    (n-type) or acceptors (p-type). Reference: Brooks-Herring; Chattopadhyay &
    Queisser, Rev. Mod. Phys. 53, 745 (1981). `energy` in J, scalar or array."""
    lambda_ac = cf.electron_mfp
    N_I = cf.doping_concentration

    if not N_I:
        # No ionized impurities: crystal MFP, with optional legacy power law.
        p = cf.electron_mfp_energy_exponent
        if p == 0:
            return lambda_ac
        return lambda_ac * (energy / (k * cf.temp)) ** p

    # Brooks-Herring ionized-impurity MFP (Matthiessen-combined with the acoustic MFP):
    energy = np.maximum(energy, 1e-4 * elementary_charge)   # guard against E -> 0
    m_star = _carrier_dos_mass(material)
    eps = material.dielectric_constant * epsilon_0
    T = cf.temp
    # Debye screening from the free carriers (n = N_I, full activation): 1/L_D^2 = e^2 n / (eps kB T)
    inv_screen_length_sq = elementary_charge**2 * N_I / (eps * k * T)
    b = 8.0 * m_star * energy / (hbar**2 * inv_screen_length_sq)     # b = (2k/beta)^2
    screening = np.maximum(np.log1p(b) - b / (1.0 + b), 1e-12)
    inv_tau_ii = (N_I * elementary_charge**4) / (16.0 * np.pi * np.sqrt(2.0 * m_star) * eps**2 * energy**1.5) * screening
    v = np.sqrt(2.0 * energy / m_star)
    lambda_ii = v / inv_tau_ii                                        # = v * tau_ii
    return 1.0 / (1.0 / lambda_ac + 1.0 / lambda_ii)


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
        With ELECTRON_MFP_ENERGY_EXPONENT != 0 the MFP becomes energy-dependent
        (see electron_mfp_at_energy); scattering stays elastic so the electron's energy,
        and hence its MFP, is fixed for the whole flight.
        """
        mfp = electron_mfp(self.energy, material)
        self.time_of_internal_scattering = -np.log(random()) * (mfp / self.speed)
