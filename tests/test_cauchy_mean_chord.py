"""
Cauchy mean-chord benchmark: in a closed cavity with fully diffuse (Lambert) walls and no
internal scattering, the mean free path between wall collisions is exactly 4V/S, for any
shape, convex or not (the radiance in the cavity is isotropic, and the collision rate is
the flux onto the walls, S*v/4 per unit density, against the V*v of path per unit time).

The particle is advanced analytically straight to the next wall (the free flight is not
what is under test, and stepping it would take minutes), and then handed to the real
scattering.surface_scattering() for the bounce. The flight lengths are measured here, so
the test checks the angular distributions of all walls TOGETHER against an exact number.
The pre-issue-#31 code gave a mean chord off by several percent.
"""

import random
from math import pi

import numpy as np
import pytest

import freepaths.scattering as SC
import freepaths.scatterers as H
from freepaths.scattering_types import ScatteringTypes

from surface_cases import EPS, Probe, make_cf, unit


N_CHORDS = 100000
TOLERANCE = 0.01     # relative; the statistical error at N_CHORDS is ~0.3%


def distance_to_box(p, d, size):
    """Distance from p (inside) along d to the walls of the box |x|<W/2, 0<y<L, |z|<T/2"""
    W, L, T = size
    t = np.inf
    for value, low, high, component in ((p[0], -W/2, W/2, d[0]), (p[1], 0.0, L, d[1]), (p[2], -T/2, T/2, d[2])):
        if component > 0:
            t = min(t, (high - value) / component)
        elif component < 0:
            t = min(t, (low - value) / component)
    return t


def distance_to_cylinder(p, d, center, radius):
    """Distance from p (outside) along d to a vertical cylinder, or inf if it is missed"""
    ox, oy = p[0] - center[0], p[1] - center[1]
    a = d[0]**2 + d[1]**2
    if a == 0:
        return np.inf
    b = 2 * (ox * d[0] + oy * d[1])
    c = ox**2 + oy**2 - radius**2
    disc = b**2 - 4 * a * c
    if disc < 0:
        return np.inf
    t = (-b - np.sqrt(disc)) / (2 * a)
    return t if t > 0 else np.inf


def mean_chord(size, hole=None, seed=31):
    """Mean distance between consecutive diffuse wall hits, in a box with an optional
    full-depth circular hole (center, radius)"""
    W, L, T = size
    cf = make_cf(width=W, length=L, thickness=T,
                 holes=[H.CircularHole(x=hole[0][0], y=hole[0][1], diameter=2 * hole[1])] if hole else None)
    SC.cf = cf
    random.seed(seed)
    rng = np.random.default_rng(seed)

    def restart():
        """Start in the middle in a random direction; this first flight is not a chord"""
        return np.array([0.0, L / 4, 0.0]), unit(rng.normal(size=3))

    p, d = restart()
    is_chord = False
    chords, restarts = [], 0
    while len(chords) < N_CHORDS:
        t = distance_to_box(p, d, size)
        if hole:
            t = min(t, distance_to_cylinder(p, d, hole[0], hole[1]))
        if is_chord:
            chords.append(t)
        hit = p + t * d
        pt = Probe(hit - EPS * d, d)
        st = ScatteringTypes()
        SC.surface_scattering(pt, st, None)
        if not (st.walls or st.top_bottom or st.holes):
            # A hit within ~EPS of an edge, where the side-wall functions reject the new
            # direction (no_new_scattering) and register nothing. Rare; start over:
            restarts += 1
            p, d = restart()
            is_chord = False
            continue
        p = np.array([pt.x, pt.y, pt.z])
        d = pt.direction()
        is_chord = True
    assert restarts < 1e-3 * N_CHORDS, f"{restarts} unregistered wall hits"
    return np.mean(chords)


def test_mean_chord_in_a_box():
    W, L, T = 300e-9, 400e-9, 100e-9
    expected = 4 * W * L * T / (2 * (W * L + W * T + L * T))
    measured = mean_chord((W, L, T))
    assert measured == pytest.approx(expected, rel=TOLERANCE), \
        f"mean chord {measured*1e9:.2f} nm, Cauchy 4V/S = {expected*1e9:.2f} nm"


def test_mean_chord_in_a_box_with_a_hole():
    W, L, T = 300e-9, 400e-9, 100e-9
    center, radius = (0.0, 200e-9), 80e-9
    volume = (W * L - pi * radius**2) * T
    surface = 2 * (W * L - pi * radius**2) + 2 * (W + L) * T + 2 * pi * radius * T
    expected = 4 * volume / surface
    measured = mean_chord((W, L, T), hole=(center, radius))
    assert measured == pytest.approx(expected, rel=TOLERANCE), \
        f"mean chord {measured*1e9:.2f} nm, Cauchy 4V/S = {expected*1e9:.2f} nm"
