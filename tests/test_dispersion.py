"""
Sanity of the tabulated phonon dispersions (materials.py stores ordinary frequency f in Hz,
not omega, in material.dispersion[:, branch + 1]; a 2 pi slip there silently scales every
velocity and kappa):
  - sound velocities at q -> 0 match the values each dispersion is built on,
  - zone-edge frequencies match the measured ones,
  - transverse branches flatten at the zone edge (v_g -> 0),
  - wavevector() never maps a finite frequency to k = 0, and inverts the tabulated points,
  - flexural ZA of graphite is quadratic (v -> 0 at q -> 0).
"""

import numpy as np
import pytest

from freepaths.materials import Si, SiGe, SiC, Diamond, Graphite


def sound_velocity(material, branch):
    d = material.dispersion
    return 2 * np.pi * (d[1, branch + 1] - d[0, branch + 1]) / (d[1, 0] - d[0, 0])


def edge_velocity(material, branch):
    d = material.dispersion
    return 2 * np.pi * abs(d[-1, branch + 1] - d[-2, branch + 1]) / (d[-1, 0] - d[-2, 0])


# (material, LA sound velocity, TA sound velocity [m/s], LA and TA zone-edge frequencies [THz]):
REFERENCE = [
    (Si, 8600, 6800, 11.3, 4.5),         # Hopkins et al. APL 95 (2009) fit; TA(X) = 4.5 THz measured
    (Diamond, 17500, 12800, 35.9, 24.2), # LA(X) ~ 35.5, TA(X) = 24.2 THz measured
    (SiGe, 7800, 6000, 13.2, 4.5),
]


@pytest.mark.parametrize("material_class, v_la, v_ta, f_la, f_ta", REFERENCE, ids=[r[0].__name__ for r in REFERENCE])
def test_sound_velocities_and_zone_edge(material_class, v_la, v_ta, f_la, f_ta):
    material = material_class(300)
    assert sound_velocity(material, 0) == pytest.approx(v_la, rel=0.03)
    assert sound_velocity(material, 1) == pytest.approx(v_ta, rel=0.03)
    assert sound_velocity(material, 2) == pytest.approx(v_ta, rel=0.03)
    assert material.dispersion[:, 1].max() / 1e12 == pytest.approx(f_la, rel=0.03)
    assert material.dispersion[:, 2].max() / 1e12 == pytest.approx(f_ta, rel=0.03)


@pytest.mark.parametrize("material_class", [Si, SiGe, SiC, Diamond])
def test_transverse_branches_flatten_at_the_zone_edge(material_class):
    material = material_class(300)
    for branch in (1, 2):
        assert edge_velocity(material, branch) < 0.1 * sound_velocity(material, branch)


def test_graphite_flexural_branch_is_quadratic():
    graphite = Graphite(300)
    za = graphite.dispersion_branch_names.index("ZA")
    assert sound_velocity(graphite, za) < 0.01 * sound_velocity(graphite, graphite.dispersion_branch_names.index("TA"))
    d = graphite.dispersion
    small = slice(1, 20)
    exponent = np.polyfit(np.log(d[small, 0]), np.log(d[small, za + 1]), 1)[0]
    assert exponent == pytest.approx(2.0, abs=0.1), f"ZA goes as q^{exponent:.2f} near Gamma"


@pytest.mark.parametrize("material_class", [Si, SiGe, SiC, Diamond, Graphite])
def test_wavevector_lookup(material_class):
    material = material_class(300)
    d = material.dispersion
    for branch in range(3):
        f_max = d[:, branch + 1].max()
        for f in np.linspace(f_max * 1e-4, f_max, 200):
            k = material.wavevector(branch, f)
            assert k > 0, f"f = {f:.3e} Hz maps onto k = 0"
            assert material.group_velocity(branch, f) > 0
        # On a monotonic branch, the tabulated points must be found exactly:
        if np.all(np.diff(d[:, branch + 1]) > 0):
            for i in range(1, len(d), 97):
                assert material.wavevector(branch, d[i, branch + 1]) == d[i, 0]
