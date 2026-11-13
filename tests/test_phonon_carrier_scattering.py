"""
Phonon-free-carrier (deformation-potential) elastic scattering, added to Si.phonon_scattering_rates
for DOPING_CONCENTRATION > 0: Asheghi, Kurabayashi, Kasnavi, Goodson, J. Appl. Phys. 91, 5079 (2002),
Eq. 11 - the degenerate ("metallic state") simplification of Ziman's free-carrier rate:

    1/tau(w) = (m* * E_D)^2 * w / (2*pi*rho*hbar^3*v),   q = w/v < 2*kF

E_D is the DILATATION deformation potential, interpolated in log10(n) over Asheghi's fitted values
(their Tables I and III) - NOT the shear potential Ju (9 eV electrons / 4.5 eV holes), which belongs
to their bound-carrier equations for the sub-Mott regime.

Checked against an independent reconstruction from scipy constants (never importing the
implementation's own prefactor), plus a bulk-RTA sanity band against the first-principles
kappa(n) of Liao et al., PRL 114, 115901 (2015), Fig. 4.
"""

import numpy as np
import pytest
from scipy.constants import electron_volt, electron_mass, hbar, k as k_B, pi

from freepaths.materials import Si, create_material


DENSITY_SI = 2330.0
SOUND_VELOCITY_SI = 6000.0
IMPURITY_COEFF = 1.8516e-45

# Asheghi Tables I + III, as (log10 n[cm^-3], E_D[eV]) and the paired fitted mass:
ASHEGHI = {
    True:  dict(mass=0.90 * electron_mass,
                log_n=[16.875, 17.398, 17.672, 18.000, 19.301, 20.230],
                e_d=[0.00, 0.05, 0.22, 0.50, 1.33, 2.33]),
    False: dict(mass=0.58 * electron_mass,
                log_n=[17.000, 18.000, 19.000],
                e_d=[0.16, 4.00, 5.00]),
}


def reference_rate(omega, doping, is_electron):
    """Independent reconstruction of Asheghi Eq. 11, not sharing code with materials.py."""
    if doping <= 0:
        return 0.0
    p = ASHEGHI[is_electron]
    e_d = np.interp(np.log10(doping * 1e-6), p["log_n"], p["e_d"]) * electron_volt
    prefactor = (p["mass"] * e_d) ** 2 / (2 * pi * DENSITY_SI * hbar ** 3 * SOUND_VELOCITY_SI)
    cutoff = 2 * (3 * pi ** 2 * doping) ** (1 / 3) * SOUND_VELOCITY_SI
    return prefactor * omega if omega < cutoff else 0.0


@pytest.fixture(scope="module")
def omegas():
    return 2 * pi * np.array([0.1e12, 0.5e12, 1e12, 3e12, 8e12])


def _carrier_only(material, omega):
    """phonon_scattering_rates returns (inelastic, elastic); strip the omega^4 impurity term."""
    return material.phonon_scattering_rates(omega)[1] - IMPURITY_COEFF * omega ** 4


def _inelastic_rate(material, omega):
    return material.phonon_scattering_rates(omega)[0]


class TestZeroDopingIsExactlyTheOldBehaviour:

    def test_undoped_matches_pre_existing_impurity_rate_exactly(self, omegas):
        si = Si(300, doping_concentration=0.0)
        for omega in omegas:
            assert si.phonon_scattering_rates(omega)[1] == pytest.approx(IMPURITY_COEFF * omega ** 4, rel=1e-12)

    def test_default_constructor_argument_is_zero_doping(self, omegas):
        si_default, si_explicit = Si(300), Si(300, doping_concentration=0.0)
        for omega in omegas:
            assert si_default.phonon_scattering_rates(omega)[1] == si_explicit.phonon_scattering_rates(omega)[1]


class TestMatchesTheClosedForm:

    @pytest.mark.parametrize("is_electron", [True, False])
    @pytest.mark.parametrize("doping", [1e24, 1e25, 2e25, 1e26])
    def test_matches_independent_reconstruction(self, omegas, is_electron, doping):
        si = Si(300, doping_concentration=doping, is_carrier_electron=is_electron)
        for omega in omegas:
            assert _carrier_only(si, omega) == pytest.approx(reference_rate(omega, doping, is_electron), rel=1e-9)

    def test_rate_is_temperature_independent(self):
        """Unlike Umklapp (exp(-theta/T)), the degenerate Ziman rate carries no explicit T."""
        omega = 2 * pi * 1e12
        rates = [_carrier_only(Si(T, doping_concentration=2e25, is_carrier_electron=False), omega)
                 for T in (10.0, 100.0, 300.0)]
        assert rates[0] == pytest.approx(rates[1], rel=1e-12) == pytest.approx(rates[2], rel=1e-12)

    def test_scales_as_E_D_squared(self):
        """The only doping dependence in the prefactor is through E_D, quadratically."""
        omega = 2 * pi * 0.05e12  # below the cutoff at both dopings
        si_lo = Si(300, doping_concentration=1e24, is_carrier_electron=False)   # 1e18 -> E_D = 4 eV
        si_hi = Si(300, doping_concentration=1e25, is_carrier_electron=False)   # 1e19 -> E_D = 5 eV
        assert _carrier_only(si_hi, omega) / _carrier_only(si_lo, omega) == pytest.approx((5.0 / 4.0) ** 2, rel=1e-6)

    def test_holes_scatter_more_strongly_than_electrons(self):
        """Liao et al. 2015's headline asymmetry: holes are the more efficient phonon scatterers
        in Si at equal carrier concentration. Asheghi's fitted parameters reproduce it."""
        omega = 2 * pi * 0.05e12
        si_p = Si(300, doping_concentration=1e25, is_carrier_electron=False)
        si_n = Si(300, doping_concentration=1e25, is_carrier_electron=True)
        assert _carrier_only(si_p, omega) > _carrier_only(si_n, omega)


class TestDeformationPotentialTable:

    def test_interpolates_between_tabulated_points(self):
        """Half-way in log10(n) between the 1e18 (4 eV) and 1e19 (5 eV) hole points."""
        omega = 2 * pi * 0.05e12
        si_mid = Si(300, doping_concentration=10 ** 24.5, is_carrier_electron=False)
        expected_e_d = 4.5  # linear in log10(n)
        si_ref = Si(300, doping_concentration=1e25, is_carrier_electron=False)  # E_D = 5 eV
        assert _carrier_only(si_mid, omega) / _carrier_only(si_ref, omega) == pytest.approx((expected_e_d / 5.0) ** 2, rel=1e-6)

    def test_clamps_above_the_tabulated_range(self):
        """np.interp clamps rather than extrapolating - Asheghi's fit must not be run off its end."""
        omega = 2 * pi * 0.05e12
        si_top = Si(300, doping_concentration=1e25, is_carrier_electron=False)      # 1e19, top of table
        si_beyond = Si(300, doping_concentration=1e27, is_carrier_electron=False)   # 1e21, far beyond
        assert _carrier_only(si_beyond, omega) == pytest.approx(_carrier_only(si_top, omega), rel=1e-9)

    def test_electron_channel_switches_off_below_the_mott_transition(self):
        """E_D -> 0 at the bottom of the electron table: carriers are bound to dopants, so the
        free-carrier channel does not exist (that regime is Asheghi's Eqs. 8/13, not modelled)."""
        omega = 2 * pi * 0.05e12
        si = Si(300, doping_concentration=1e22, is_carrier_electron=True)  # 1e16 cm^-3, below table
        assert _carrier_only(si, omega) == pytest.approx(0.0, abs=1e-30)


class TestHardCutoff:

    def test_rate_is_exactly_zero_above_cutoff(self):
        si = Si(300, doping_concentration=2e25, is_carrier_electron=False)
        cutoff = si._phonon_carrier_cutoff_omega
        assert _carrier_only(si, cutoff * 1.5) == pytest.approx(0.0, abs=1e-30)
        assert _carrier_only(si, cutoff * 0.5) > 0.0

    def test_cutoff_is_2kF_times_sound_velocity(self):
        doping = 2e25
        si = Si(300, doping_concentration=doping, is_carrier_electron=False)
        expected = 2 * (3 * pi ** 2 * doping) ** (1 / 3) * SOUND_VELOCITY_SI
        assert si._phonon_carrier_cutoff_omega == pytest.approx(expected, rel=1e-9)

    def test_higher_doping_pushes_the_cutoff_up(self):
        """kF ~ n^(1/3): more doping extends scattering to higher omega, not just a bigger rate."""
        si_lo = Si(300, doping_concentration=1e25, is_carrier_electron=False)
        si_hi = Si(300, doping_concentration=8e25, is_carrier_electron=False)
        assert si_hi._phonon_carrier_cutoff_omega == pytest.approx(2 * si_lo._phonon_carrier_cutoff_omega, rel=1e-9)


class TestBulkKappaAgainstLiaoFirstPrinciples:
    """Sanity band on the integrated effect, against Liao et al., PRL 114, 115901 (2015), Fig. 4.
    Ziman's classical formula is expected to OVERpredict the reduction (Liao call it 'apparently
    insufficient for a modern understanding'), so these are generous bands, not tight fits -
    they exist to catch an implementation that is orders of magnitude off, which is exactly the
    failure mode of a mistranscribed prefactor."""

    @staticmethod
    def kappa_rta(material, temp=300.0):
        disp = material.dispersion
        k = disp[:, 0]
        total = 0.0
        for branch in (1, 2, 3):
            f = disp[:, branch]
            valid = f > 1e9
            kk, omega = k[valid], 2 * pi * f[valid]
            vg = np.gradient(omega, kk)
            x = np.clip(hbar * omega / (k_B * temp), 1e-8, 700)
            c = k_B * x ** 2 * np.exp(x) / np.expm1(x) ** 2
            total += np.trapezoid(kk ** 2 * c * vg ** 2 * material.phonon_relaxation_time(omega), kk) / (2 * pi ** 2) / 3
        return total

    def test_undoped_reproduces_the_documented_rta_value(self):
        assert self.kappa_rta(Si(300.0, doping_concentration=0.0)) == pytest.approx(157.7, rel=0.02)

    @pytest.mark.parametrize("doping_cm3, lo, hi", [
        (1e18, 0.00, 0.25),   # Liao: negligible below 1e18
        (1e19, 0.05, 0.40),   # Liao p-type ~13%
        (1e21, 0.30, 0.85),   # Liao p-type ~47%
    ])
    def test_hole_reduction_is_in_a_physically_sane_band(self, doping_cm3, lo, hi):
        k0 = self.kappa_rta(Si(300.0, doping_concentration=0.0))
        kd = self.kappa_rta(Si(300.0, doping_concentration=doping_cm3 * 1e6, is_carrier_electron=False))
        reduction = 1 - kd / k0
        assert lo <= reduction <= hi, f"{reduction:.1%} outside [{lo:.0%}, {hi:.0%}] at {doping_cm3:.0e} cm^-3"

    def test_holes_reduce_kappa_more_than_electrons(self):
        k_p = self.kappa_rta(Si(300.0, doping_concentration=1e26, is_carrier_electron=False))
        k_n = self.kappa_rta(Si(300.0, doping_concentration=1e26, is_carrier_electron=True))
        assert k_p < k_n


class TestInelasticChannelIsUntouched:

    def test_umklapp_rate_independent_of_doping(self, omegas):
        si_undoped = Si(300, doping_concentration=0.0)
        si_doped = Si(300, doping_concentration=2e25, is_carrier_electron=False)
        for omega in omegas:
            assert _inelastic_rate(si_doped, omega) == pytest.approx(_inelastic_rate(si_undoped, omega), rel=1e-12)


class TestFactoryWiring:

    def test_create_material_forwards_doping_to_si(self, omegas):
        si = create_material("Si", 300, doping_concentration=2e25, is_carrier_electron=False)
        for omega in omegas:
            assert _carrier_only(si, omega) == pytest.approx(reference_rate(omega, 2e25, False), rel=1e-9)

    def test_create_material_ignores_doping_for_materials_without_the_parameter(self):
        assert create_material("SiC", 300, doping_concentration=2e25, is_carrier_electron=False).name == "SiC"
