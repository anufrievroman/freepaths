"""
Diffuse (Lambert) scattering must be a cosine-law hemisphere about the true surface normal,
for every surface type the code can scatter on.

Each sample places an electron (always diffuse, since specularity() returns 0 for electrons)
just in front of a random point P on one surface, aimed at it, and calls the real dispatcher:
scattering.surface_scattering(), reinitialization() for the hot sides, or Source.generate_angles().
The outgoing direction d is compared with the normal n at P, which surface_cases.py computes
from the shape (it does not reuse any angle from the code under test). For a Lambert lobe:
  - no direction goes into the wall (mu = d.n >= 0),
  - mu has CDF mu^2, so <mu> = 2/3,
  - the tangential components average to zero (the lobe is not tilted),
  - the azimuth around n is uniform.

The pre-issue-#31 formulas fail the mu test on every surface with KS p-values ~1e-60 at 20k samples.
"""

import random
from math import pi

import numpy as np
import pytest
from scipy.stats import kstest

import freepaths.scattering as SC
from freepaths.scattering_primitives import lambert_angles
from freepaths.scattering_types import ScatteringTypes
from freepaths.sources import Source, Distributions

from surface_cases import EPS, Probe, angles_to_direction, build_cases, incoming_direction, tangent_frame


N = 20000            # samples per surface
P_MIN = 1e-4         # KS p-value threshold (seeds are fixed, so the outcome is deterministic)

CASES = build_cases()


def scatter_samples(cf, sampler, kind, rng):
    """Scatter N particles on the surface; returns {branch: list of (d, n)}"""
    SC.cf = cf
    samples = {}
    for _ in range(N):
        P, n = sampler(rng)
        d_in = incoming_direction(n, rng)
        pt = Probe(P - EPS * d_in, d_in)
        st = ScatteringTypes()
        if kind == "reinit":
            SC.reinitialization(pt, st)
            event = st.hot_side
        else:
            SC.surface_scattering(pt, st, None)
            event = st.holes or st.walls or st.top_bottom or st.pillars or st.interfaces
        assert event is not None, f"no scattering registered at P={P}, n={n}"

        branch, n_out = "", n
        if kind == "interface":
            transmitted = st.interfaces_transmission is not None
            branch = "transmitted" if transmitted else "reflected"
            n_out = -n if transmitted else n   # a transmitted lobe points onward, into the far side
        samples.setdefault(branch, []).append((pt.direction(), n_out))
    return samples


def source_samples(angle):
    source = Source(angle_distribution=Distributions.LAMBERT, angle=angle)
    n = np.array([np.sin(angle), np.cos(angle), 0.0])
    return {"": [(angles_to_direction(*source.generate_angles()), n) for _ in range(N)]}


def assert_lambert(samples, label):
    mu, t1s, t2s = [], [], []
    for d, n in samples:
        t1, t2 = tangent_frame(n)
        mu.append(d @ n)
        t1s.append(d @ t1)
        t2s.append(d @ t2)
    mu, t1s, t2s = map(np.array, (mu, t1s, t2s))
    count = len(mu)

    into_wall = np.sum(mu < -1e-9)
    assert into_wall == 0, f"{label}: {into_wall} of {count} directions point into the surface"

    p_mu = kstest(np.clip(mu, 0, 1)**2, "uniform").pvalue
    assert p_mu > P_MIN, f"{label}: cos(angle to normal) is not Lambert, <mu> = {mu.mean():.4f} (expected 2/3), KS p = {p_mu:.2g}"

    tangential_se = np.sqrt(np.mean(t1s**2 + t2s**2) / 2 / count)
    tilt = np.degrees(np.arctan2(np.hypot(t1s.mean(), t2s.mean()), mu.mean()))
    assert abs(t1s.mean()) < 5 * tangential_se and abs(t2s.mean()) < 5 * tangential_se, \
        f"{label}: lobe is tilted by {tilt:.2f} deg from the surface normal"

    p_azimuth = kstest((np.arctan2(t2s, t1s) + pi) / (2 * pi), "uniform").pvalue
    assert p_azimuth > P_MIN, f"{label}: azimuth around the normal is not uniform, KS p = {p_azimuth:.2g}"


@pytest.mark.parametrize("name, cf, sampler, kind", CASES, ids=[c[0] for c in CASES])
def test_diffuse_scattering_is_lambert_about_the_surface_normal(name, cf, sampler, kind):
    random.seed(31)
    rng = np.random.default_rng(31)
    samples = source_samples(sampler) if kind == "source" else scatter_samples(cf, sampler, kind, rng)
    if kind == "interface":
        assert set(samples) == {"reflected", "transmitted"}
    for branch, branch_samples in samples.items():
        assert_lambert(branch_samples, f"{name} {branch}".strip())


@pytest.mark.parametrize("normal_theta, normal_phi", [
    (0.0, 0.0), (pi / 2, 0.0), (-2.5, 0.0), (0.0, pi / 2), (1.0, -pi / 2), (0.8, 0.6), (-2.0, -1.1)])
def test_lambert_angles_about_arbitrary_normals(normal_theta, normal_phi):
    """lambert_angles() alone, including normals out of the xy plane (tilted pillar walls)"""
    random.seed(31)
    n = angles_to_direction(normal_theta, normal_phi)
    samples = [(angles_to_direction(*lambert_angles(normal_theta, normal_phi)), n) for _ in range(N)]
    assert_lambert(samples, f"normal ({normal_theta}, {normal_phi})")
