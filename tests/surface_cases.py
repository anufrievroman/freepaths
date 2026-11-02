"""
Shared geometry for the surface-scattering tests: a catalogue of every surface the
code can scatter on, each with a sampler that returns a point P on the surface and the
normal n at P (pointing into the region the particle comes from). The normals are
computed here from the shape itself, never from the code under test.
"""

from math import pi
from types import SimpleNamespace

import numpy as np

from freepaths.options import ParticleType
import freepaths.scatterers as H


EPS = 1e-11          # particle starts EPS before the surface and would end EPS behind it [m]
SPEED = 1e3
T, W, L = 100e-9, 4e-6, 4e-6   # thickness, width, length of the test domain [m]
X0, Y0, SX, SY = 0.0, 2e-6, 600e-9, 500e-9   # position and size of the test holes [m]


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


class Probe:
    """Minimal particle. An electron always scatters diffusely (specularity() returns 0);
    a phonon on surfaces of zero roughness always scatters specularly (it returns 1)."""
    k = 1e9
    f = 1e12
    wavelength = 2 * pi / 1e9

    def __init__(self, position, direction, particle_type=ParticleType.ELECTRON):
        self.type = particle_type
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


def build_cases():
    """
    List of (name, cf, sampler, kind):
      kind "surface"   -- scattered through scattering.surface_scattering()
      kind "reinit"    -- hot side, through scattering.reinitialization() (always diffuse)
      kind "interface" -- VerticalPlane, transmission forced to 0.5 (both branches)
      kind "source"    -- Lambert particle source; sampler is the source angle
    """
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
            # P must lie on the outer boundary, and away from the cusps where two circles meet
            # (there the normal is undefined and the nearest point can switch within EPS):
            distances, _ = line.tree.query(P, k=2)
            if abs(distances[0] - line_radius) < 1e-12 and distances[1] > line_radius + 1e-9:
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

    # Vertical interface, with the transmission forced to 0.5 so both branches are well sampled,
    # and equal group velocities on both sides, so a specular transmission goes straight through:
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
