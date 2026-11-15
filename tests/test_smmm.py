"""
SMMM interface transmission (Ran & Cao, PRB 110, 024302 (2024)), interface_smmm.py:
  - alpha_spec is the oblique-incidence acoustic mismatch 4 W_i W_j / (W_i + W_j)^2 with
    W = Z / cos(theta), i.e. crossed cosines (an uncrossed denominator was a real bug, fixed 2026-09-25);
  - it is reciprocal, vanishes past the critical angle, and is 1 for identical media;
  - alpha_diff is the DMM split by k^2, complementary between the two directions;
  - alpha_total_2T stays in [0, 1], reduces to the specular and diffuse limits, vanishes
    under total internal reflection, and its diffuse plateau gives R >= 3/4.
"""

import math

import numpy as np
import pytest

from freepaths.interface_smmm import alpha_spec, alpha_diff, alpha_total_2T, get_k_from_omega
from freepaths.materials import Si, SiGe


rng = np.random.default_rng(31)


def snell_angle(theta_i, vg_i, vg_j):
    """sin(theta_i) / vg_i = sin(theta_j) / vg_j"""
    return math.asin(vg_j / vg_i * math.sin(theta_i))


def test_alpha_spec_at_normal_incidence():
    for rho_i, rho_j, vg_i, vg_j in [(2330, 3008, 8400, 7700), (1000, 4000, 3000, 6000)]:
        z_i, z_j = rho_i * vg_i, rho_j * vg_j
        assert alpha_spec(0.0, vg_i, vg_j, rho_i, rho_j) == pytest.approx(4 * z_i * z_j / (z_i + z_j)**2, rel=1e-12)


def test_alpha_spec_is_the_oblique_acoustic_mismatch():
    for _ in range(2000):
        rho_i, rho_j = rng.uniform(1000, 6000, 2)
        vg_i, vg_j = rng.uniform(2000, 9000, 2)
        theta_i = rng.uniform(0, math.pi / 2 * 0.99)
        if vg_j / vg_i * math.sin(theta_i) >= 1:
            continue
        theta_j = snell_angle(theta_i, vg_i, vg_j)
        w_i = rho_i * vg_i / math.cos(theta_i)
        w_j = rho_j * vg_j / math.cos(theta_j)
        assert alpha_spec(theta_i, vg_i, vg_j, rho_i, rho_j) == pytest.approx(4 * w_i * w_j / (w_i + w_j)**2, rel=1e-10)


def test_alpha_spec_is_reciprocal_and_bounded():
    for _ in range(2000):
        rho_i, rho_j = rng.uniform(1000, 6000, 2)
        vg_i, vg_j = rng.uniform(2000, 9000, 2)
        theta_i = rng.uniform(0, math.pi / 2 * 0.99)
        if vg_j / vg_i * math.sin(theta_i) >= 1:
            assert alpha_spec(theta_i, vg_i, vg_j, rho_i, rho_j) == 0.0     # total internal reflection
            continue
        forward = alpha_spec(theta_i, vg_i, vg_j, rho_i, rho_j)
        backward = alpha_spec(snell_angle(theta_i, vg_i, vg_j), vg_j, vg_i, rho_j, rho_i)
        assert 0.0 <= forward <= 1.0
        assert forward == pytest.approx(backward, rel=1e-9)


def test_alpha_spec_identical_media_transmit_fully():
    for theta in (0.0, 0.5, 1.2):
        assert alpha_spec(theta, 8000, 8000, 2330, 2330) == pytest.approx(1.0, rel=1e-12)


def test_alpha_diff_is_a_complementary_k_squared_split():
    for _ in range(200):
        k_i, k_j = rng.uniform(1e8, 1e10, 2)
        p = rng.uniform(0, 0.99)
        forward = alpha_diff(k_i, k_j, p, p)
        assert forward == pytest.approx(k_j**2 / (k_i**2 + k_j**2), rel=1e-12)
        assert forward + alpha_diff(k_j, k_i, p, p) == pytest.approx(1.0, rel=1e-12)
    assert alpha_diff(2e9, 2e9, 0.3, 0.3) == pytest.approx(0.5)


@pytest.fixture(scope="module")
def si_sige():
    return Si(300), SiGe(300)


def total(theta, f, branch, roughness, si, sige):
    omega = 2 * math.pi * f
    vg_i = si.group_velocity(branch, f)
    vg_j = sige.group_velocity(branch, f)
    return alpha_total_2T(theta, vg_i, vg_j, si.density, sige.density, omega, omega, si, sige, branch, roughness)


def test_alpha_total_is_bounded(si_sige):
    si, sige = si_sige
    for _ in range(3000):
        value = total(rng.uniform(0, math.pi / 2 * 0.999), rng.uniform(0.05e12, 4e12), int(rng.integers(3)),
                      rng.uniform(0, 5e-9), si, sige)
        assert 0.0 <= value <= 1.0


def test_alpha_total_specular_limit(si_sige):
    """Smooth interface: T = a_spec(Si->SiGe) * a_spec(SiGe->Si)"""
    si, sige = si_sige
    for f, branch, theta in [(1e12, 0, 0.0), (1e12, 0, 0.6), (2e12, 1, 0.3)]:
        vg_i, vg_j = si.group_velocity(branch, f), sige.group_velocity(branch, f)
        theta_t = snell_angle(theta, vg_i, vg_j)
        expected = (alpha_spec(theta, vg_i, vg_j, si.density, sige.density)
                    * alpha_spec(theta_t, vg_j, vg_i, sige.density, si.density))
        assert total(theta, f, branch, 0.0, si, sige) == pytest.approx(expected, rel=1e-12)


def test_alpha_total_diffuse_plateau(si_sige):
    """Very rough interface: T = x(1 - x) with x = k_j^2 / (k_i^2 + k_j^2), so R = 1 - T >= 3/4"""
    si, sige = si_sige
    for f, branch in [(2e12, 0), (3e12, 1), (6e12, 0)]:
        omega = 2 * math.pi * f
        k_i = get_k_from_omega(omega, si, branch)
        k_j = get_k_from_omega(omega, sige, branch)
        x = k_j**2 / (k_i**2 + k_j**2)
        value = total(0.2, f, branch, 50e-9, si, sige)
        assert value == pytest.approx(x * (1 - x), rel=1e-9)
        assert 1 - value >= 0.75 - 1e-12


def test_total_reflection_only_on_the_way_into_the_faster_medium():
    """Snell's law written out with numbers, not through snell_angle(): from 3000 to 6000 m/s the
    critical angle is arcsin(3000/6000) = 30 deg; from 6000 to 3000 m/s there is none (the phonon
    bends towards the normal). The inverted law, sin(theta_j) = (vg_i/vg_j) sin(theta_i), has it
    the other way round."""
    slow, fast, rho = 3000.0, 6000.0, 2500.0
    assert alpha_spec(math.radians(29.9), slow, fast, rho, rho) > 0
    for degrees in (30.1, 45, 60, 89):
        assert alpha_spec(math.radians(degrees), slow, fast, rho, rho) == 0.0
    for degrees in (0, 30, 45, 60, 80, 89):
        assert alpha_spec(math.radians(degrees), fast, slow, rho, rho) > 0
    # Refraction at 30 deg from fast to slow: sin(theta_j) = 0.25
    theta_j = math.asin(0.25)
    w_i, w_j = rho * fast / math.cos(math.radians(30)), rho * slow / math.cos(theta_j)
    assert alpha_spec(math.radians(30), fast, slow, rho, rho) == pytest.approx(4 * w_i * w_j / (w_i + w_j)**2, rel=1e-12)


def test_a_slower_layer_transmits_at_every_angle(si_sige):
    """SiGe is the slower medium for these modes, so a phonon from Si is never totally reflected
    by a SiGe layer, however grazing: it bends towards the normal on the way in, and is back at
    its incidence angle on the way out."""
    si, sige = si_sige
    for f, branch in [(1e12, 0), (3e12, 0), (1e12, 1), (3e12, 1)]:
        assert sige.group_velocity(branch, f) < si.group_velocity(branch, f)
        for theta in np.linspace(0, math.pi / 2 * 0.999, 40):
            assert total(theta, f, branch, 0.5e-9, si, sige) > 0


def test_a_faster_layer_reflects_totally_past_the_critical_angle(si_sige):
    """The same pair the other way round: a Si layer in SiGe, critical angle arcsin(v_SiGe / v_Si)"""
    si, sige = si_sige
    f, branch = 1e12, 0
    omega = 2 * math.pi * f
    v_sige, v_si = sige.group_velocity(branch, f), si.group_velocity(branch, f)
    critical = math.asin(v_sige / v_si)

    def through_si_layer(theta):
        return alpha_total_2T(theta, v_sige, v_si, sige.density, si.density, omega, omega, sige, si, branch, 0.5e-9)
    assert through_si_layer(critical * 0.98) > 0
    for theta in np.linspace(critical * 1.001, math.pi / 2 * 0.999, 20):
        assert through_si_layer(theta) == 0.0


def test_no_transmission_above_the_top_of_the_branch_in_the_layer(si_sige):
    """A mode with no counterpart in the layer (frequency above the top of the branch there) cannot enter it"""
    si, sige = si_sige
    for branch in range(3):
        top_si, top_sige = si.dispersion[:, branch + 1].max(), sige.dispersion[:, branch + 1].max()
        if top_si > top_sige:
            f = (top_si + top_sige) / 2
            assert total(0.3, f, branch, 0.5e-9, si, sige) == 0.0
        else:
            f = (top_si + top_sige) / 2
            omega = 2 * math.pi * f
            assert alpha_total_2T(0.3, sige.group_velocity(branch, f), si.group_velocity(branch, f), sige.density,
                                  si.density, omega, omega, sige, si, branch, 0.5e-9) == 0.0
