"""
Diffuse (Lambert) scattering must be a cosine-law hemisphere about the true surface normal,
for every surface type the code can scatter on.

Each sample places an electron (always diffuse, since specularity() returns 0 for electrons)
just in front of a random point P on one surface, aimed at it, and calls the real dispatcher:
scattering.surface_scattering(), reinitialization() for the hot sides, or Source.generate_angles().
The outgoing direction d is compared with the normal n at P, which the test computes on its own
from the shape (it does not reuse any angle from the code under test). For a Lambert lobe:
  - no direction goes into the wall (mu = d.n >= 0),
  - mu has CDF mu^2, so <mu> = 2/3,
  - the tangential components average to zero (the lobe is not tilted),
  - the azimuth around n is uniform.

The pre-issue-#31 formulas fail the mu test on every surface with KS p-values ~1e-130 at 20k samples.
"""

import random
from math import pi
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.stats import kstest

from freepaths.options import ParticleType
import freepaths.scattering as SC
import freepaths.scatterers as H
from freepaths.scattering_primitives import lambert_angles
from freepaths.scattering_types import ScatteringTypes
from freepaths.sources import Source, Distributions


N = 20000            # samples per surface
P_MIN = 1e-4         # KS p-value threshold (seeds are fixed, so the outcome is deterministic)
EPS = 1e-11          # particle starts EPS before the surface and would end EPS behind it [m]
SPEED = 1e3
T, W, L = 100e-9, 4e-6, 4e-6   # thickness, width, length of the test domain [m]
X0, Y0, SX, SY = 0.0, 2e-6, 600e-9, 500e-9   # position and size of the test holes [m]


# ---------------------------------------------------------------- helpers

def make_cf(**kwargs):
    """Minimal config for the scattering dispatcher, with all roughness 0 and no structures"""
    cf = SimpleNamespace(
        timestep=2 * EPS / SPEED, thickness=T, width=W, length=L, depletion_width=0.0,
        hole_roughness=0.0, side_wall_roughness=0.0, top_roughness=0.0, pillar_roughness=0.0,
        interface_roughness=0.0, include_right_sidewall=True, include_left_sidewall=True,
        include_top_sidewall=True, include_bottom_sidewall=True, holes=None, pillars=None,
        interfaces=None, temp=300, is_carrier_electron=True,
        hot_side_position_bottom=False, hot_side_position_top=False,
        hot_side_position_right=False, hot_side_position_left=False)
    cf.__dict__.update(kwargs)
    return cf


class Electron:
    """Minimal particle: always scatters diffusely"""
    type = ParticleType.ELECTRON
    k = 1e9
    f = 1e12

    def __init__(self, position, direction):
        self.x, self.y, self.z = position
        self.theta = np.arctan2(direction[0], direction[1])
        self.phi = np.arcsin(np.clip(direction[2], -1, 1))
        self.speed = SPEED

    def speed_in(self, material):
        return SPEED

    def correct_angle(self):
        if abs(self.theta) > pi:
            self.theta -= np.sign(self.theta) * 2 * pi

    def direction(self):
        return angles_to_direction(self.theta, self.phi)


def angles_to_direction(theta, phi):
    """Direction convention of move.py"""
    return np.array([np.sin(theta) * np.cos(phi), np.cos(theta) * np.cos(phi), np.sin(phi)])


def unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def tangent_frame(n):
    t1 = unit(np.cross([0, 0, 1], n)) if abs(n[2]) < 0.99 else np.array([1.0, 0, 0])
    return t1, np.cross(n, t1)


def incoming_direction(n, rng):
    """Random direction towards the surface, avoiding near-grazing incidence that makes
    the face choice ambiguous within EPS of an edge"""
    t1, t2 = tangent_frame(n)
    mu = np.sqrt(rng.uniform(0.01, 1.0))
    psi = 2 * pi * rng.random()
    s = np.sqrt(1 - mu**2)
    return -(mu * n + s * np.cos(psi) * t1 + s * np.sin(psi) * t2)


def z_random(rng, margin=0.05):
    return (rng.random() - 0.5) * T * (1 - 2 * margin)


def segment_sampler(a, b, n):
    """Points on segment a-b in the xy plane (away from its ends) with the in-plane normal n"""
    a, b, n = np.asarray(a, float), np.asarray(b, float), unit([n[0], n[1], 0.0])

    def sample(rng):
        xy = a + rng.uniform(0.05, 0.95) * (b - a)
        return np.array([xy[0], xy[1], z_random(rng)]), n
    return sample


def triangle_edges(vertices):
    """Edges of a triangle as (a, b, normal pointing away from the triangle)"""
    center = np.mean(vertices, axis=0)
    edges = []
    for i in range(3):
        a, b = np.asarray(vertices[i], float), np.asarray(vertices[(i + 1) % 3], float)
        n = unit([(b - a)[1], -(b - a)[0]])
        if np.dot(n, (a + b) / 2 - center) < 0:
            n = -n
        edges.append((a, b, n))
    return edges


def edge_label(n):
    if abs(n[1]) > 0.999:
        return "base"
    if abs(n[0]) > 0.999:
        return "vertical cut"
    return "slant"


# ---------------------------------------------------------------- surfaces

def build_cases():
    """(name, cf, sampler, kind): sampler(rng) returns a point P on the surface and the
    normal at P pointing into the region the particle comes from"""
    cases = []

    def across_width(rng):
        return W * (rng.random() - 0.5) * 0.8

    def along_length(rng):
        return L * rng.uniform(0.1, 0.9)

    # Domain walls:
    cases += [
        ("side wall right", make_cf(), lambda r: (np.array([W/2, along_length(r), z_random(r)]), np.array([-1.0, 0, 0])), "surface"),
        ("side wall left", make_cf(), lambda r: (np.array([-W/2, along_length(r), z_random(r)]), np.array([1.0, 0, 0])), "surface"),
        ("end wall top", make_cf(), lambda r: (np.array([across_width(r), L, z_random(r)]), np.array([0, -1.0, 0])), "surface"),
        ("end wall bottom", make_cf(), lambda r: (np.array([across_width(r), 0.0, z_random(r)]), np.array([0, 1.0, 0])), "surface"),
        ("ceiling", make_cf(), lambda r: (np.array([across_width(r), along_length(r), T/2]), np.array([0, 0, -1.0])), "surface"),
        ("floor", make_cf(), lambda r: (np.array([across_width(r), along_length(r), -T/2]), np.array([0, 0, 1.0])), "surface"),
    ]

    # Hot-side re-emission:
    hot_sides = {
        "bottom": (lambda r: [across_width(r), 0.0, z_random(r)], [0, 1, 0]),
        "top": (lambda r: [across_width(r), L, z_random(r)], [0, -1, 0]),
        "right": (lambda r: [W/2, along_length(r), z_random(r)], [-1, 0, 0]),
        "left": (lambda r: [-W/2, along_length(r), z_random(r)], [1, 0, 0]),
    }
    for side, (point, n) in hot_sides.items():
        cases.append((f"hot side {side}", make_cf(**{f"hot_side_position_{side}": True}),
                      lambda r, point=point, n=n: (np.array(point(r), float), np.array(n, float)), "reinit"))

    # Circular hole:
    radius = 200e-9

    def circle(r):
        a = 2 * pi * r.random()
        n = np.array([np.cos(a), np.sin(a), 0.0])
        return np.array([X0 + radius * n[0], Y0 + radius * n[1], z_random(r)]), n
    cases.append(("circular hole", make_cf(holes=[H.CircularHole(x=X0, y=Y0, diameter=2 * radius)]), circle, "surface"))

    # Rectangular hole, each face, and the floor of a partial-depth hole hit from below:
    hx, hy = SX / 2, SY / 2
    faces = {
        "left face": ((X0 - hx, Y0 - hy), (X0 - hx, Y0 + hy), (-1, 0)),
        "right face": ((X0 + hx, Y0 - hy), (X0 + hx, Y0 + hy), (1, 0)),
        "top face": ((X0 - hx, Y0 + hy), (X0 + hx, Y0 + hy), (0, 1)),
        "bottom face": ((X0 - hx, Y0 - hy), (X0 + hx, Y0 - hy), (0, -1)),
    }
    for face, (a, b, n) in faces.items():
        cases.append((f"rectangular hole {face}", make_cf(holes=[H.RectangularHole(x=X0, y=Y0, size_x=SX, size_y=SY)]),
                      segment_sampler(a, b, n), "surface"))
    depth = 40e-9
    cases.append(("partial rectangular hole floor",
                  make_cf(holes=[H.RectangularHole(x=X0, y=Y0, size_x=SX, size_y=SY, depth=depth)]),
                  lambda r: (np.array([X0 + SX * (r.random() - 0.5) * 0.9, Y0 + SY * (r.random() - 0.5) * 0.9, T/2 - depth]),
                             np.array([0, 0, -1.0])), "surface"))

    # Triangular holes (up = apex at +y), full and both halves:
    tri_up = [(X0 - SX/2, Y0 - SY/2), (X0 + SX/2, Y0 - SY/2), (X0, Y0 + SY/2)]
    tri_down = [(X0 - SX/2, Y0 + SY/2), (X0 + SX/2, Y0 + SY/2), (X0, Y0 - SY/2)]
    for name, cls, vertices in [("triangular up hole", H.TriangularUpHole, tri_up),
                                ("triangular down hole", H.TriangularDownHole, tri_down)]:
        for a, b, n in triangle_edges(vertices):
            label = edge_label(n) if edge_label(n) == "base" else ("right slant" if n[0] > 0 else "left slant")
            cases.append((f"{name} {label}", make_cf(holes=[cls(x=X0, y=Y0, size_x=SX, size_y=SY)]),
                          segment_sampler(a, b, n), "surface"))
    for name, cls, full in [("half triangular down hole", H.TriangularDownHalfHole, tri_down),
                            ("half triangular up hole", H.TriangularUpHalfHole, tri_up)]:
        for is_right in (True, False):
            vertices = [(X0, full[0][1]), full[1] if is_right else full[0], full[2]]
            for a, b, n in triangle_edges(vertices):
                cases.append((f"{name} ({'right' if is_right else 'left'}) {edge_label(n)}",
                              make_cf(holes=[cls(x=X0, y=Y0, size_x=SX, size_y=SY, is_right_half=is_right)]),
                              segment_sampler(a, b, n), "surface"))

    # Point-line hole (union of circles along a line); the normal is radial from the nearest point:
    line_radius = 40e-9
    line = H.PointLineHole(x=X0, y=Y0, points=[(i * 10e-9, 0.0) for i in range(-30, 31)], thickness=2 * line_radius)

    def point_line(r):
        while True:
            center = line.points[r.integers(len(line.points))]
            a = 2 * pi * r.random()
            n = np.array([np.cos(a), np.sin(a)])
            P = center + line_radius * n
            if abs(line.tree.query(P)[0] - line_radius) < 1e-12:   # P lies on the outer boundary
                return np.array([P[0], P[1], z_random(r)]), np.array([n[0], n[1], 0.0])
    cases.append(("point-line hole", make_cf(holes=[line]), point_line, "surface"))

    # Parabolic walls: top is y = -x^2/(4f) + tip (hole above), bottom is y = x^2/(4f) + tip (hole below):
    focus = 1e-6

    def parabola(is_top, tip):
        def sample(r):
            x = across_width(r)
            x = x if abs(x) > 1e-9 else 1e-9
            if is_top:
                return np.array([x, -x**2 / (4 * focus) + tip, z_random(r)]), unit([-x / (2 * focus), -1.0, 0.0])
            return np.array([x, x**2 / (4 * focus) + tip, z_random(r)]), unit([-x / (2 * focus), 1.0, 0.0])
        return sample
    cases.append(("parabola top wall", make_cf(holes=[H.ParabolaTop(tip=3e-6, focus=focus)]), parabola(True, 3e-6), "surface"))
    cases.append(("parabola bottom wall", make_cf(holes=[H.ParabolaBottom(tip=1e-6, focus=focus)]), parabola(False, 1e-6), "surface"))

    # Pillar wall, hit from inside the pillar; the normal points to the pillar axis:
    pillar_radius, pillar_height = 150e-9, 300e-9

    def pillar(r):
        a = 2 * pi * r.random()
        radial = np.array([np.cos(a), np.sin(a), 0.0])
        z = T/2 + pillar_height * r.uniform(0.1, 0.9)
        return np.array([X0 + pillar_radius * radial[0], Y0 + pillar_radius * radial[1], z]), -radial
    cases.append(("pillar wall", make_cf(pillars=[H.CircularPillar(x=X0, y=Y0, diameter=2 * pillar_radius, height=pillar_height)]),
                  pillar, "surface"))

    # Vertical interface, with the transmission forced to 0.5 so both branches are well sampled:
    interface = H.VerticalPlane(position_x=0.0, outer_material="Si", inner_material="SiGe")
    interface.assign_materials(SimpleNamespace(temp=300))
    interface.vg_i_to_vg_j = 1.0
    interface.transmission_probability = lambda cf, pt: 0.5

    def interface_point(r):
        side = 1.0 if r.random() < 0.5 else -1.0
        return np.array([0.0, along_length(r), z_random(r)]), np.array([side, 0, 0])
    cases.append(("vertical interface", make_cf(interfaces=[interface]), interface_point, "interface"))

    # Lambert particle sources:
    for angle in (0.0, 0.7, -2.0):
        cases.append((f"source angle {angle}", None, angle, "source"))

    return cases


CASES = build_cases()


# ---------------------------------------------------------------- sampling and checks

def scatter_samples(cf, sampler, kind, rng):
    """Scatter N particles on the surface; returns {branch: list of (d, n)}"""
    SC.cf = cf
    samples = {}
    for _ in range(N):
        P, n = sampler(rng)
        d_in = incoming_direction(n, rng)
        pt = Electron(P - EPS * d_in, d_in)
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
