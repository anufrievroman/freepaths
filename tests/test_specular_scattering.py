"""
Specular reflection must obey the mirror law about the true surface normal, d' = d - 2(d.n)n,
on every surface that can reflect specularly.

A phonon on a surface of zero roughness has specularity exactly 1, so the specular branch
always runs. As in the Lambert test, the particle is aimed at a random point P of the surface
through the real dispatcher, and the outgoing direction is compared with the mirror image of
the incoming one about the normal at P, computed from the geometry in surface_cases.py.

For the interface, a specular transmission with equal group velocities on both sides
(vg_i/vg_j = 1, no refraction) must leave the direction unchanged.
"""

import random

import numpy as np
import pytest

import freepaths.scattering as SC
from freepaths.options import ParticleType
from freepaths.scattering_types import ScatteringTypes, Scattering

from surface_cases import EPS, Probe, build_cases, incoming_direction


N = 2000
TOLERANCE = 2e-3     # the curved surfaces take their normal from the would-be position EPS behind P

CASES = [c for c in build_cases() if c[3] in ("surface", "interface")]


@pytest.mark.parametrize("name, cf, sampler, kind", CASES, ids=[c[0] for c in CASES])
def test_specular_reflection_obeys_the_mirror_law(name, cf, sampler, kind):
    random.seed(31)
    rng = np.random.default_rng(31)
    SC.cf = cf
    branches = set()
    for _ in range(N):
        P, n = sampler(rng)
        d_in = incoming_direction(n, rng)
        pt = Probe(P - EPS * d_in, d_in, ParticleType.PHONON)
        st = ScatteringTypes()
        SC.surface_scattering(pt, st, None)
        event = st.holes or st.walls or st.top_bottom or st.pillars or st.interfaces
        assert event is Scattering.SPECULAR, f"{name}: expected a specular event at P={P}, got {event}"

        d_out = pt.direction()
        if kind == "interface" and st.interfaces_transmission is not None:
            branches.add("transmitted")
            expected, what = d_in, "transmission without refraction"
        else:
            branches.add("reflected")
            expected, what = d_in - 2 * (d_in @ n) * n, "mirror reflection"
        error = np.linalg.norm(d_out - expected)
        assert error < TOLERANCE, (f"{name}: {what} is wrong at P={P}, n={n}: "
                                   f"in {d_in}, out {d_out}, expected {expected}")
        assert abs(np.linalg.norm(d_out) - 1) < 1e-9

    if kind == "interface":
        assert branches == {"reflected", "transmitted"}
