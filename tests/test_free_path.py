"""
Free-path bookkeeping of the real particle loop (run_particle), which the MFP-sampling kappa
is built on. Ground truth is counted independently of the Flight: every particle.move() is
recorded, so the distance actually travelled is (number of moves) x (step length).
  - Specular events (walls, interfaces) must NOT end a free path: a particle that only meets
    specular surfaces has ONE free path, equal to the whole distance travelled.
  - Diffuse events must end it: the number of saved free paths is (diffuse events + 1), and
    the free paths add up to the whole distance.
  - The internal-scattering clock must advance with every step the particle moves.
"""

import random

import numpy as np
import pytest

import freepaths.scatterers as H
from freepaths.options import ParticleType

from particle_loop import SPEED, TIMESTEP, LoopParticle, configure, run
from surface_cases import unit


STEP = SPEED * TIMESTEP


def launch_directions(n, rng):
    """Directions with a positive y component, so the particle heads for the cold side"""
    directions = []
    while len(directions) < n:
        d = unit(rng.normal(size=3))
        if d[1] > 0.3:
            directions.append(d)
    return directions


def assert_one_free_path(flight, particle, label):
    moves = len(particle.positions)
    assert particle.y > 0, f"{label}: the particle should have reached the cold side"
    assert flight._mfp_count == 1, f"{label}: {flight._mfp_count} free paths, expected 1 (no diffuse event happened)"
    assert flight._mfp_sum == pytest.approx(moves * STEP, rel=1e-9), \
        f"{label}: free path {flight._mfp_sum*1e9:.1f} nm, distance travelled {moves*STEP*1e9:.1f} nm"
    assert flight.time_since_previous_scattering == pytest.approx(moves * TIMESTEP, rel=1e-9), \
        f"{label}: scattering clock {flight.time_since_previous_scattering:.3e} s, time travelled {moves*TIMESTEP:.3e} s"


def test_specular_walls_do_not_end_the_free_path(monkeypatch):
    configure(monkeypatch)
    random.seed(31)
    rng = np.random.default_rng(31)
    for d in launch_directions(30, rng):
        particle = LoopParticle((0.0, 1e-9, 0.0), d, ParticleType.PHONON)
        flight = run(particle, monkeypatch)
        assert_one_free_path(flight, particle, "specular walls")


def test_specular_interface_transmission_does_not_end_the_free_path(monkeypatch):
    """Interfaces every 14 nm (the Layers geometry), always transmitted specularly without refraction"""
    interfaces = []
    for x in np.arange(-86e-9, 86e-9, 14e-9):
        plane = H.VerticalPlane(position_x=float(x), outer_material="Si", inner_material="SiGe")
        plane.transmission_probability = lambda cf, pt: 1.0
        plane.vg_i_to_vg_j = 1.0
        interfaces.append(plane)
    configure(monkeypatch, interfaces=interfaces)
    random.seed(31)
    rng = np.random.default_rng(31)
    for d in launch_directions(30, rng):
        particle = LoopParticle((1e-9, 1e-9, 0.0), d, ParticleType.PHONON)
        flight = run(particle, monkeypatch)
        assert_one_free_path(flight, particle, "specular interface transmission")


def test_diffuse_events_end_the_free_path(monkeypatch):
    configure(monkeypatch, length=400e-9)
    random.seed(31)
    rng = np.random.default_rng(31)
    total_diffuse = 0
    for d in launch_directions(30, rng):
        particle = LoopParticle((0.0, 1e-9, 0.0), d, ParticleType.ELECTRON)
        flight, diffuse = run(particle, monkeypatch, count_diffuse=True)
        total_diffuse += diffuse
        moves = len(particle.positions)
        assert particle.y > 400e-9 or particle.y < 0, "the particle should have left through the cold or hot side"
        assert flight._mfp_count == diffuse + 1, f"{flight._mfp_count} free paths for {diffuse} diffuse events"
        assert flight._mfp_sum == pytest.approx(moves * STEP, rel=1e-9), \
            f"free paths add up to {flight._mfp_sum*1e9:.1f} nm, distance travelled {moves*STEP*1e9:.1f} nm"
    assert total_diffuse > 100, "the walls should have scattered the electrons many times"
