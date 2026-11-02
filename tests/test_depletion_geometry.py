"""
Surface-depletion dead layer of every hole shape, Hole.is_inside_depleted(x, y, z, cf, W):
  - with W = 0 it is identical to is_inside,
  - with W > 0 it is a superset of is_inside,
  - it extends exactly W beyond every face: a point W - delta away from the hole along the
    face normal is in the dead layer, a point W + delta away is not.
The faces and their normals come from surface_cases.py (computed from the shapes).
"""

import numpy as np
import pytest

from surface_cases import build_cases, make_cf


DELTA = 0.2e-9
WIDTHS = [3e-9, 7e-9]
N = 400

CASES = [c for c in build_cases() if c[3] == "surface" and c[1].holes and "parabola" not in c[0]]


def outward_points(sampler, rng):
    """Points P on a hole face with the normal n pointing away from the hole"""
    for _ in range(N):
        yield sampler(rng)


@pytest.mark.parametrize("name, cf, sampler, kind", CASES, ids=[c[0] for c in CASES])
def test_dead_layer_geometry(name, cf, sampler, kind):
    hole = cf.holes[0]
    rng = np.random.default_rng(31)
    for P, n in outward_points(sampler, rng):
        # Random points around the face, inside and outside the hole:
        for s in rng.uniform(-20e-9, 20e-9, 5):
            x, y, z = P + s * n
            assert hole.is_inside_depleted(x, y, z, cf, 0.0) == hole.is_inside(x, y, z, cf)
            for width in WIDTHS:
                if hole.is_inside(x, y, z, cf):
                    assert hole.is_inside_depleted(x, y, z, cf, width), f"{name}: dead layer must contain the hole"
        # The dead layer extends exactly W along the face normal:
        for width in WIDTHS:
            if "point-line" in name:
                # The union of circles: only where the offset point is still closest to its own circle
                distances, _ = hole.tree.query((P + (width + DELTA) * n)[:2], k=1)
                if abs(distances - (hole.thickness / 2 + width + DELTA)) > 1e-12:
                    continue
            assert hole.is_inside_depleted(*(P + (width - DELTA) * n), cf, width), \
                f"{name}: point {(width - DELTA)*1e9:.1f} nm off the face should be in the {width*1e9:.0f} nm dead layer"
            assert not hole.is_inside_depleted(*(P + (width + DELTA) * n), cf, width), \
                f"{name}: point {(width + DELTA)*1e9:.1f} nm off the face should be outside the {width*1e9:.0f} nm dead layer"
