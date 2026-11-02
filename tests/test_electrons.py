"""
Electron physics and post-processing, checked against closed forms:
  - parabolic-band speed, Brooks-Herring ionized-impurity MFP (Matthiessen with the acoustic MFP),
    surface-depletion width;
  - BTE conductivity, Seebeck coefficient and electronic thermal conductivity for a constant MFP,
    where Xi(E) is linear in E and the integrals reduce to Fermi-Dirac integrals
    G_n = int x^n f dx: sigma = e^2 A kT G0, S = -/+ (k/e) (2 G1/G0 - eta),
    kappa_e = A k^3 T^2 (3 G2 - 4 eta G1 + eta^2 G0) - S^2 sigma T;
  - the MC-to-BTE calibration: a diffusive 1/ToF ~ v reproduces the BTE exactly, a fixed
    mapping constant scales sigma but not S;
  - the kappa-weighted specularity accumulators (SpecularityData) against a direct per-phonon sum.
"""

from types import SimpleNamespace

import numpy as np
import pytest
from scipy.constants import electron_volt, elementary_charge as e, epsilon_0, hbar, k as k_B
from scipy.integrate import quad

from freepaths.config import cf
from freepaths.electron import Electron, depletion_width, electron_mfp
from freepaths.materials import Si
from freepaths.post_computations import ElectronPostComputation
from freepaths.data import SpecularityData
from freepaths.flight import Flight, BOUNDARY_CATEGORIES
from freepaths.scattering_types import ScatteringTypes, Scattering


@pytest.fixture
def config(monkeypatch):
    """Electron settings on the shared config, restored after each test"""
    settings = dict(temp=300.0, media="Si", media_fermi_level=0.0, is_carrier_electron=True,
                    electron_mfp=15e-9, doping_concentration=0, electron_mfp_energy_exponent=0,
                    depletion_doping=0, surface_potential=0.0, mean_mapping_constant=None,
                    energy_step=1e-3, fermi_level_lower_bound=-0.2, fermi_level_upper_bound=0.2)
    for name, value in settings.items():
        monkeypatch.setattr(cf, name, value, raising=False)
    return cf


@pytest.fixture(scope="module")
def silicon():
    return Si(300)


# ---------------------------------------------------------------- particle-level physics

@pytest.mark.parametrize("is_electron", [True, False])
def test_parabolic_band_speed(silicon, is_electron):
    mass = silicon.effective_electron_dos_mass if is_electron else silicon.effective_hole_dos_mass
    for energy in (0.01, 0.1, 0.3):
        particle = SimpleNamespace(energy=energy * electron_volt, is_electron_carrier=is_electron)
        assert Electron.speed_in(particle, silicon) == pytest.approx(np.sqrt(2 * energy * electron_volt / mass), rel=1e-12)


def test_mfp_without_doping_is_the_crystal_mfp(config, silicon):
    energies = np.linspace(0.005, 0.3, 20) * electron_volt
    assert np.all(electron_mfp(energies, silicon) == config.electron_mfp)


def brooks_herring_rate(energy, density, material, temp):
    """1/tau of Brooks-Herring screened-Coulomb scattering, as written in Lundstrom,
    Fundamentals of Carrier Transport, Eq. (2.86)"""
    m = material.effective_electron_dos_mass
    eps = material.dielectric_constant * epsilon_0
    debye_length_sq = eps * k_B * temp / (e**2 * density)
    gamma_sq = 8 * m * energy * debye_length_sq / hbar**2
    return (density * e**4 / (16 * np.sqrt(2 * m) * np.pi * eps**2) * energy**-1.5
            * (np.log(1 + gamma_sq) - gamma_sq / (1 + gamma_sq)))


def test_brooks_herring_mfp(config, monkeypatch, silicon):
    density = 1e25    # 1e19 cm^-3
    monkeypatch.setattr(cf, "doping_concentration", density)
    energies = np.linspace(0.005, 0.3, 30) * electron_volt
    mfp = electron_mfp(energies, silicon)
    speed = np.sqrt(2 * energies / silicon.effective_electron_dos_mass)
    impurity_mfp = speed / brooks_herring_rate(energies, density, silicon, cf.temp)
    assert np.allclose(1 / mfp, 1 / cf.electron_mfp + 1 / impurity_mfp, rtol=1e-10)   # Matthiessen
    assert np.all(mfp < cf.electron_mfp)
    assert np.all(np.diff(mfp) > 0), "impurities must scatter low-energy carriers most"

    # More dopants mean a shorter MFP only where screening is weak (b >> 1). At low energy the
    # Debye screening, which also tightens with N, wins: 1/tau ~ E^(1/2) / N for b << 1, so at
    # 1e20 cm^-3 a 5 meV carrier has a LONGER MFP than at 1e19. That is Brooks-Herring itself
    # (nondegenerate Debye screening), not a coding error:
    monkeypatch.setattr(cf, "doping_concentration", 1e26)
    weak_screening = energies > 0.05 * electron_volt
    assert np.all(electron_mfp(energies, silicon)[weak_screening] < mfp[weak_screening]), \
        "more dopants must mean a shorter MFP in the weak-screening range"


def test_depletion_width(config, monkeypatch, silicon):
    assert depletion_width(silicon) == 0.0                       # no surface potential
    monkeypatch.setattr(cf, "surface_potential", 0.4)
    assert depletion_width(silicon) == 0.0                       # no doping
    monkeypatch.setattr(cf, "doping_concentration", 1e25)
    expected = np.sqrt(2 * 11.7 * epsilon_0 * 0.4 / (e * 1e25))
    assert depletion_width(silicon) == pytest.approx(expected, rel=1e-12)
    assert depletion_width(silicon) == pytest.approx(7.19e-9, rel=1e-3)
    monkeypatch.setattr(cf, "depletion_doping", 4e25)            # overrides DOPING_CONCENTRATION
    assert depletion_width(silicon) == pytest.approx(expected / 2, rel=1e-12)


# ---------------------------------------------------------------- post-processing

def run_post_processing(time_of_flight, n_per_energy=3, n_lost=50, energy_top=0.8):
    """ElectronPostComputation on synthetic data: every energy of the grid (1 meV steps up to
    energy_top eV) gets the same ToF for each of its electrons, plus some electrons that never
    arrived (ToF = 0)"""
    energies = np.arange(1, round(energy_top * 1000) + 1) * 1e-3 * electron_volt
    initial = np.repeat(energies, n_per_energy)
    times = np.repeat(time_of_flight(energies), n_per_energy)
    initial = np.concatenate((initial, energies[:n_lost]))
    times = np.concatenate((times, np.zeros(n_lost)))
    post = ElectronPostComputation(SimpleNamespace(initial_energies=initial, travel_times=times))
    post.compute()
    return post


def fermi_integrals(eta):
    """G_n = int_0^inf x^n / (1 + exp(x - eta)) dx for n = 0, 1, 2"""
    upper = max(eta, 0) + 60
    return [quad(lambda x: x**n / (1 + np.exp(x - eta)), 0, upper)[0] for n in range(3)]


@pytest.mark.parametrize("is_electron", [True, False])
def test_bte_conductivity_and_seebeck_match_closed_forms(config, monkeypatch, is_electron):
    monkeypatch.setattr(cf, "is_carrier_electron", is_electron)
    post = run_post_processing(lambda E: np.ones_like(E))
    material = post.material
    m_dos = material.effective_electron_dos_mass if is_electron else material.effective_hole_dos_mass
    m_cond = material.effective_electron_mass if is_electron else material.effective_hole_mass
    # Xi(E) = (mfp/3) * sqrt(2E/m_cond) * g(E) = A * E:
    A = cf.electron_mfp / 3 * np.sqrt(2 / m_cond) * (2 * m_dos / hbar**2)**1.5 / (2 * np.pi**2)
    kT = k_B * cf.temp
    sign = -1 if is_electron else 1
    for row in range(0, len(post.fermi_levels), 40):
        eta = post.fermi_levels[row] / kT
        g0, g1, g2 = fermi_integrals(eta)
        sigma = e**2 * A * kT * g0
        seebeck = sign * (k_B / e) * (2 * g1 / g0 - eta)
        kappa_e = A * k_B**3 * cf.temp**2 * (3 * g2 - 4 * eta * g1 + eta**2 * g0) - seebeck**2 * sigma * cf.temp
        assert post.bte_conductivity[row, 1] == pytest.approx(sigma, rel=1e-2)
        assert post.bte_seebeck[row, 1] == pytest.approx(seebeck, rel=1e-2)
        assert np.sign(post.bte_seebeck[row, 1]) == sign
        assert post.bte_thermal_conductivity[row, 1] == pytest.approx(kappa_e, rel=1e-2)
        assert post.bte_thermal_conductivity[row, 1] > 0


def test_lorenz_number_limits(config, monkeypatch):
    """L = kappa_e / (sigma T) = (k/e)^2 (3 G2/G0 - (2 G1/G0)^2) for Xi ~ E: it tends to 2 (k/e)^2 in
    the nondegenerate limit (G_n = n! e^eta) and to the Wiedemann-Franz pi^2/3 (k/e)^2 deep in the
    degenerate one"""
    monkeypatch.setattr(cf, "fermi_level_upper_bound", 0.6)
    post = run_post_processing(lambda E: np.ones_like(E), energy_top=1.6)
    lorenz = post.bte_thermal_conductivity[:, 1] / (post.bte_conductivity[:, 1] * cf.temp)
    unit = (k_B / e)**2
    eta = post.fermi_levels / (k_B * cf.temp)
    assert lorenz[np.argmin(abs(eta + 6))] == pytest.approx(2 * unit, rel=1e-2)
    assert lorenz[np.argmin(abs(eta - 20))] == pytest.approx(np.pi**2 / 3 * unit, rel=2e-2)


def test_diffusive_monte_carlo_reproduces_the_bte(config, silicon):
    """1/ToF ~ v(E) is the diffusive bulk result; the auto-calibrated MC must then equal the BTE"""
    speed = lambda E: np.sqrt(2 * E / silicon.effective_electron_mass)
    post = run_post_processing(lambda E: 1e-9 / speed(E))
    assert np.allclose(post.mc_conductivity[:, 1], post.bte_conductivity[:, 1], rtol=1e-9)
    assert np.allclose(post.mc_seebeck[:, 1], post.bte_seebeck[:, 1], rtol=1e-9)
    assert np.allclose(post.mc_thermal_conductivity[:, 1], post.bte_thermal_conductivity[:, 1], rtol=1e-9)


def test_fixed_mapping_constant_scales_sigma_but_not_seebeck(config, monkeypatch, silicon):
    speed = lambda E: np.sqrt(2 * E / silicon.effective_electron_mass)
    tof = lambda E: 1e-9 / speed(E) * (1 + 0.3 * E / electron_volt)   # not exactly diffusive
    auto = run_post_processing(tof)
    monkeypatch.setattr(cf, "mean_mapping_constant", 2 * auto.mean_mapping_constant)
    fixed = run_post_processing(tof)
    assert np.allclose(fixed.mc_conductivity[:, 1], 2 * auto.mc_conductivity[:, 1], rtol=1e-12)
    assert np.allclose(fixed.mc_seebeck[:, 1], auto.mc_seebeck[:, 1], rtol=1e-12)


# ---------------------------------------------------------------- specularity accumulators

def test_boundary_events_are_tallied_by_category():
    flight = SimpleNamespace(boundary_events=[0] * (2 * len(BOUNDARY_CATEGORIES)))
    st = ScatteringTypes()
    st.walls, st.top_bottom, st.holes, st.interfaces = (Scattering.SPECULAR, Scattering.DIFFUSE,
                                                         Scattering.DIFFUSE, Scattering.SPECULAR)
    st.interfaces_transmission = Scattering.SPECULAR    # a subset of interfaces, must not count twice
    Flight.record_boundary_events(flight, st)
    assert flight.boundary_events == [1, 0, 0, 1, 0, 1, 0, 0, 1, 0]


def test_specularity_accumulators_match_the_per_phonon_sum():
    rng = np.random.default_rng(31)
    n_categories = len(BOUNDARY_CATEGORIES)
    flights = []
    for _ in range(3000):
        events = rng.integers(0, 20, size=2 * n_categories) * (rng.random(2 * n_categories) < 0.6)
        flights.append(SimpleNamespace(boundary_events=list(events), thermal_conductivity=rng.exponential(1.0)))

    # Three workers, merged like main_mfp_sampling does:
    merged = SpecularityData()
    for worker_flights in np.array_split(np.arange(len(flights)), 3):
        worker = SpecularityData()
        for i in worker_flights:
            worker.add_particle(flights[i])
        merged.read_data(worker.dump_data())

    events = np.array([f.boundary_events for f in flights], float).reshape(len(flights), n_categories, 2)
    kappa = np.array([f.thermal_conductivity for f in flights])
    specular, diffuse = events[:, :, 0], events[:, :, 1]
    total = specular + diffuse
    assert np.allclose(merged.specular_events, specular.sum(0))
    assert np.allclose(merged.diffuse_events, diffuse.sum(0))
    for c in range(n_categories):
        hit = total[:, c] > 0
        weighted = np.sum(kappa[hit] * specular[hit, c] / total[hit, c]) / np.sum(kappa[hit])
        assert merged.kappa_weighted_specular[c] / merged.kappa_of_scattered[c] == pytest.approx(weighted, rel=1e-10)
