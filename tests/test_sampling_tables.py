"""
Dispersion-based phonon sampling (Peraud & Hadjiconstantinou, PRB 84, 205331 (2011)):
  - hot-wall emission is weighted by D(w) C(w) v(w) per dispersion bin,
  - rethermalization at inelastic events by D(w) C(w) / tau_inelastic(w),
so that both populations are C(w)-shaped once the residence time is accounted for,
  - draw_from_table() samples those tables and returns k and v of exactly the drawn bin,
  - an internal event is attributed to each channel in proportion to its rate.
The heat capacity here is computed by differentiating the Bose-Einstein occupation
numerically, independently of the closed form used in materials.py.
"""

import random
from math import pi
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.constants import hbar, k as k_B
from scipy.stats import chisquare

from freepaths.materials import Si, SiGe, SiC
from freepaths.phonon import Phonon


@pytest.fixture(scope="module", params=[10, 300])
def silicon(request):
    return Si(request.param)


def bins(material, branch):
    """Bin-midpoint frequency, k^2 dk weight and bin-averaged group velocity, from the raw dispersion"""
    k_grid = material.dispersion[:, 0]
    f_grid = material.dispersion[:, branch + 1]
    f_mid = (f_grid[1:] + f_grid[:-1]) / 2
    keep = f_mid > 0
    dos = (((k_grid[1:] + k_grid[:-1]) / 2)**2 * np.diff(k_grid))[keep]
    velocity = (2 * pi * np.abs(np.diff(f_grid)) / np.diff(k_grid))[keep]
    return f_mid[keep], dos, velocity


def heat_capacity(f, temp):
    """Mode heat capacity hbar*w * dn/dT by central finite difference of the Bose-Einstein occupation"""
    def occupation(t):
        return 1 / np.expm1(hbar * 2 * pi * f / (k_B * t))
    dt = temp * 1e-4
    return hbar * 2 * pi * f * (occupation(temp + dt) - occupation(temp - dt)) / (2 * dt)


def bin_probabilities(cumulative):
    return np.diff(np.concatenate(([0.0], cumulative)))


def test_cumulative_tables_are_valid(silicon):
    for table in (silicon.emission_frequency_probabilities, silicon.scattering_frequency_probabilities):
        for cumulative in table:
            assert np.all(np.diff(cumulative) >= 0)
            assert cumulative[-1] == pytest.approx(1.0, abs=1e-12)
    for cumulative in (silicon.emission_branch_probabilities, silicon.scattering_branch_probabilities):
        assert len(cumulative) == 3 and np.all(np.diff(cumulative) >= 0)
        assert cumulative[-1] == pytest.approx(1.0, abs=1e-12)


def test_tables_hold_the_dispersion_bins(silicon):
    for branch in range(3):
        f, _, velocity = bins(silicon, branch)
        assert np.allclose(silicon.frequencies_table[branch], f, rtol=1e-12)
        assert np.allclose(silicon.group_velocity_table[branch], velocity, rtol=1e-12)
        assert np.all(silicon.wavevector_table[branch] > 0)


@pytest.mark.parametrize("table", ["emission", "scattering"])
def test_bin_weights(silicon, table):
    """Joint (branch, bin) probability must be proportional to D*C*v (emission) or D*C/tau_inelastic (scattering)"""
    branch_cumulative = getattr(silicon, f"{table}_branch_probabilities")
    branch_probability = bin_probabilities(branch_cumulative)
    ratios = []
    for branch in range(3):
        f, dos, velocity = bins(silicon, branch)
        if table == "emission":
            weight = dos * heat_capacity(f, silicon.temp) * velocity
        else:
            inelastic = np.array([silicon.phonon_scattering_rates(2 * pi * x)[0] for x in f])
            weight = dos * heat_capacity(f, silicon.temp) * inelastic
        probability = branch_probability[branch] * bin_probabilities(getattr(silicon, f"{table}_frequency_probabilities")[branch])
        significant = weight > weight.max() * 1e-6     # far tails lose relative precision in both
        ratios.append(probability[significant] / weight[significant])
    ratios = np.concatenate(ratios)
    assert np.allclose(ratios, ratios[0], rtol=1e-5), \
        f"{table} probabilities are not proportional to the expected weight (spread {ratios.std()/ratios.mean():.2e})"


def test_steady_populations_are_heat_capacity_shaped(silicon):
    """Emission / v and retherm * tau_inelastic must both give the same D*C population"""
    emission_branch = bin_probabilities(silicon.emission_branch_probabilities)
    scattering_branch = bin_probabilities(silicon.scattering_branch_probabilities)
    emitted, rethermalized = [], []
    for branch in range(3):
        f, _, velocity = bins(silicon, branch)
        inelastic = np.array([silicon.phonon_scattering_rates(2 * pi * x)[0] for x in f])
        emitted.append(emission_branch[branch] * bin_probabilities(silicon.emission_frequency_probabilities[branch]) / velocity)
        rethermalized.append(scattering_branch[branch] * bin_probabilities(silicon.scattering_frequency_probabilities[branch]) / inelastic)
    emitted, rethermalized = np.concatenate(emitted), np.concatenate(rethermalized)
    emitted, rethermalized = emitted / emitted.sum(), rethermalized / rethermalized.sum()
    significant = emitted > emitted.max() * 1e-6
    assert np.allclose(emitted[significant], rethermalized[significant], rtol=1e-5)


def test_draw_from_table_samples_the_tables():
    material = Si(300)
    random.seed(31)
    n_draws = 100000
    n_groups = 20
    counts = np.zeros((3, n_groups))
    for _ in range(n_draws):
        phonon = SimpleNamespace()
        Phonon.draw_from_table(phonon, material, material.emission_branch_probabilities,
                               material.emission_frequency_probabilities)
        branch = phonon.branch_number
        index = int(np.flatnonzero(material.frequencies_table[branch] == phonon.f)[0])
        # k and v must belong to exactly the drawn bin:
        assert phonon.k == material.wavevector_table[branch][index]
        assert phonon.speed == material.group_velocity_table[branch][index]
        counts[branch, index * n_groups // len(material.frequencies_table[branch])] += 1

    branch_probability = bin_probabilities(material.emission_branch_probabilities)
    expected = np.zeros((3, n_groups))
    for branch in range(3):
        per_bin = bin_probabilities(material.emission_frequency_probabilities[branch])
        groups = np.arange(len(per_bin)) * n_groups // len(per_bin)
        expected[branch] = branch_probability[branch] * np.bincount(groups, weights=per_bin, minlength=n_groups)
    expected *= n_draws
    keep = expected > 5
    p = chisquare(counts[keep], expected[keep] * counts[keep].sum() / expected[keep].sum()).pvalue
    assert p > 1e-4, f"draws do not follow the emission table, chi-square p = {p:.2g}"


@pytest.mark.parametrize("material_class", [Si, SiGe, SiC])
def test_relaxation_time_is_the_sum_of_the_channel_rates(material_class):
    material = material_class(300)
    for f in (0.5e12, 2e12, 6e12):
        omega = 2 * pi * f
        inelastic, elastic = material.phonon_scattering_rates(omega)
        assert inelastic >= 0 and elastic >= 0
        assert material.phonon_relaxation_time(omega) == pytest.approx(1 / (inelastic + elastic), rel=1e-12)


def test_internal_event_channel_split(monkeypatch):
    """At an internal event, the inelastic channel must be chosen with probability rate_inelastic / rate_total"""
    import freepaths.phonon as phonon_module
    monkeypatch.setattr(phonon_module.cf, "phonon_hydrodynamic", False)
    material = SiGe(300)          # SiGe: a large elastic (alloy) share
    random.seed(31)
    for f in (1e12, 5e12):
        inelastic, elastic = material.phonon_scattering_rates(2 * pi * f)
        expected = inelastic / (inelastic + elastic)
        phonon = SimpleNamespace(f=f, _grain_boundary_rate=lambda: 0.0)
        n = 40000
        hits = sum(Phonon.classify_internal_event(phonon, material) == "inelastic" for _ in range(n))
        sigma = np.sqrt(expected * (1 - expected) / n)
        assert abs(hits / n - expected) < 5 * sigma, f"f = {f:.1e}: inelastic share {hits/n:.4f}, expected {expected:.4f}"


def test_ballistic_conductance():
    """Landauer G = (1/4) sum over branches of int d^3k/(2 pi)^3 C v: for one speed v it is C v / 4"""
    import gray_material
    gray = gray_material.Gray(300)
    assert gray.ballistic_conductance() == pytest.approx(gray.dispersion_heat_capacity * gray_material.VELOCITY / 4, rel=1e-9)
    silicon = Si(300)
    total = 0.0
    for branch in range(3):
        f, dos, velocity = bins(silicon, branch)
        total += np.sum(dos / (2 * pi**2) * heat_capacity(f, silicon.temp) * velocity) / 4
    assert silicon.ballistic_conductance() == pytest.approx(total, rel=1e-6)
