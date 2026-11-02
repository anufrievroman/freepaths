"""
Containment: in the real particle loop, a particle must never leave the structure or enter
a hole. Every position after every move is checked against the geometry: inside the box
(or inside a pillar, which extends the structure above the ceiling) and outside every hole.
Only the very last position may lie beyond the cold or hot side, since that is how a
particle exits. Electrons (always diffuse) and phonons (a mix of specular and diffuse on
rough walls) are run through a membrane with holes of every shape, a pillar and an interface.
"""

import random

import numpy as np
import pytest

import freepaths.scatterers as H
from freepaths.config import cf
from freepaths.options import ParticleType

from particle_loop import LoopParticle, configure, run
from surface_cases import unit


W, L, T = 600e-9, 800e-9, 100e-9
PILLAR = dict(x=-150e-9, y=150e-9, diameter=80e-9, height=100e-9)


def structure():
    interface = H.VerticalPlane(position_x=250e-9, outer_material="Si", inner_material="SiGe")
    interface.transmission_probability = lambda cf, pt: 0.5
    interface.vg_i_to_vg_j = 1.0
    holes = [
        H.CircularHole(x=-150e-9, y=400e-9, diameter=120e-9),
        H.RectangularHole(x=100e-9, y=250e-9, size_x=120e-9, size_y=90e-9),
        H.RectangularHole(x=100e-9, y=650e-9, size_x=120e-9, size_y=90e-9, depth=40e-9),
        H.TriangularUpHole(x=-150e-9, y=650e-9, size_x=120e-9, size_y=100e-9),
        H.TriangularDownHole(x=100e-9, y=450e-9, size_x=120e-9, size_y=100e-9),
        H.TriangularDownHalfHole(x=-50e-9, y=550e-9, size_x=100e-9, size_y=80e-9, is_right_half=True),
        H.PointLineHole(x=150e-9, y=120e-9, points=[(i * 10e-9, 0.0) for i in range(-3, 4)], thickness=30e-9),
    ]
    return holes, [H.CircularPillar(**PILLAR)], [interface]


def in_pillar(x, y, z):
    return (np.hypot(x - PILLAR["x"], y - PILLAR["y"]) <= PILLAR["diameter"] / 2
            and T / 2 <= z <= T / 2 + PILLAR["height"])


def violations(positions, holes):
    """Positions outside the structure or inside a hole (the last position may be beyond y = 0 or y = L)"""
    bad = []
    for i, (x, y, z) in enumerate(positions):
        last = i == len(positions) - 1
        inside_box = abs(x) <= W / 2 and abs(z) <= T / 2 and (last or 0 <= y <= L)
        if not (inside_box or in_pillar(x, y, z)):
            bad.append((i, (x, y, z), "outside the structure"))
        elif any(hole.is_inside(x, y, z, cf) for hole in holes):
            bad.append((i, (x, y, z), "inside a hole"))
    return bad


@pytest.mark.parametrize("particle_type, roughness", [(ParticleType.ELECTRON, 0.0), (ParticleType.PHONON, 0.3e-9)],
                         ids=["electrons", "phonons"])
def test_particles_stay_inside_the_structure(monkeypatch, particle_type, roughness):
    holes, pillars, interfaces = structure()
    configure(monkeypatch, width=W, length=L, thickness=T, holes=holes, pillars=pillars, interfaces=interfaces,
              number_of_timesteps=60000, side_wall_roughness=roughness, top_roughness=roughness,
              hole_roughness=roughness, pillar_roughness=roughness, interface_roughness=roughness)
    random.seed(31)
    rng = np.random.default_rng(31)
    failures = []
    for n in range(150):
        while True:
            start = np.array([rng.uniform(-0.48, 0.48) * W, 2e-9, rng.uniform(-0.45, 0.45) * T])
            if not any(hole.is_inside(*start, cf) for hole in holes):
                break
        d = unit(rng.normal(size=3))
        d[1] = abs(d[1])
        particle = LoopParticle(start, d, particle_type)
        run(particle, monkeypatch)
        bad = violations(particle.positions, holes)
        if bad:
            failures.append((n, len(particle.positions), bad[0]))
    assert not failures, (f"{len(failures)} of 150 particles left the structure or entered a hole; "
                          f"first ones (particle, steps, (step, position, what)): {failures[:5]}")
