"""
Soffer specularity p = exp(-16 pi^2 sigma^2 cos^2(a) / lambda^2):
  - the formula, its limits, and lambda = 2 pi / k for phonons;
  - every surface must pass the angle a between the incoming direction and its TRUE normal
    (cos a = |d.n|), since an angle measured from the wrong axis silently makes grazing
    phonons diffuse or normal ones specular.
"""

import random
from math import pi
from types import SimpleNamespace

import numpy as np
import pytest

import freepaths.scattering as SC
import freepaths.scattering_primitives as primitives
import freepaths.scattering_interfaces as interfaces
import freepaths.scatterers as scatterers
from freepaths.options import ParticleType
from freepaths.phonon import Phonon
from freepaths.scattering_primitives import specularity
from freepaths.scattering_types import ScatteringTypes

from surface_cases import EPS, Probe, build_cases, incoming_direction


def phonon_with_k(k):
    return SimpleNamespace(type=ParticleType.PHONON, wavelength=2 * pi / k)


def test_phonon_wavelength_is_two_pi_over_k():
    for k in (1e7, 3.3e8, 1.1e10):
        assert Phonon.wavelength.fget(SimpleNamespace(k=k)) == pytest.approx(2 * pi / k, rel=1e-15)


@pytest.mark.parametrize("angle", [0.0, 0.3, 1.0, 1.4])
@pytest.mark.parametrize("roughness", [0.1e-9, 0.5e-9, 2e-9])
@pytest.mark.parametrize("k", [5e8, 2e9, 1e10])
def test_soffer_formula(angle, roughness, k):
    expected = np.exp(-16 * pi**2 * roughness**2 * np.cos(angle)**2 * k**2 / (4 * pi**2))
    assert specularity(angle, roughness, phonon_with_k(k)) == pytest.approx(expected, rel=1e-12)


def test_soffer_limits():
    pt = phonon_with_k(2e9)
    assert specularity(0.4, 0.0, pt) == 1.0                       # perfectly smooth surface
    assert specularity(pi / 2, 5e-9, pt) == pytest.approx(1.0)    # grazing incidence
    assert specularity(0.0, 5e-9, pt) < 1e-6                      # rough surface, short wavelength
    assert specularity(0.0, 1e-9, phonon_with_k(1e5)) == pytest.approx(1.0, abs=1e-6)   # long wavelength
    electron = SimpleNamespace(type=ParticleType.ELECTRON)
    assert specularity(pi / 2, 0.0, electron) == 0                # electrons are always diffuse


def test_soffer_is_monotonic():
    pt = phonon_with_k(2e9)
    angles = np.linspace(0, pi / 2, 50)
    p_angle = [specularity(a, 0.5e-9, pt) for a in angles]
    assert np.all(np.diff(p_angle) >= 0), "p must grow towards grazing incidence"
    roughnesses = np.linspace(0, 3e-9, 50)
    p_roughness = [specularity(0.2, s, pt) for s in roughnesses]
    assert np.all(np.diff(p_roughness) <= 0), "p must fall with roughness"
    ks = np.linspace(1e8, 1e10, 50)
    p_k = [specularity(0.2, 0.5e-9, phonon_with_k(k)) for k in ks]
    assert np.all(np.diff(p_k) <= 0), "p must fall with wavevector"


CASES = [c for c in build_cases() if c[3] in ("surface", "interface")]
N = 2000


@pytest.mark.parametrize("name, cf, sampler, kind", CASES, ids=[c[0] for c in CASES])
def test_surfaces_pass_the_angle_to_their_normal(name, cf, sampler, kind, monkeypatch):
    """The angle handed to specularity() must satisfy cos(a) = |d.n| at the hit point"""
    angles = []

    def recording_specularity(angle, roughness, particle):
        angles.append(angle)
        return 1.0

    for module in (primitives, interfaces, scatterers):
        monkeypatch.setattr(module, "specularity", recording_specularity)

    random.seed(31)
    rng = np.random.default_rng(31)
    SC.cf = cf
    worst = (0.0, None)
    for _ in range(N):
        P, n = sampler(rng)
        d_in = incoming_direction(n, rng)
        pt = Probe(P - EPS * d_in, d_in, ParticleType.PHONON)
        angles.clear()
        SC.surface_scattering(pt, ScatteringTypes(), None)
        assert len(angles) == 1, f"{name}: specularity() called {len(angles)} times"
        error = abs(abs(np.cos(angles[0])) - abs(d_in @ n))
        if error > worst[0]:
            worst = (error, (P, n, d_in, angles[0]))
    assert worst[0] < 2e-3, (f"{name}: angle passed to specularity() is not the angle to the normal; "
                             f"worst |cos a| - |d.n| = {worst[0]:.3f} at (P, n, d_in, a) = {worst[1]}")
