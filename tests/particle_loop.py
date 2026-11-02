"""
Helpers to run the real run_particle() loop on a test particle in a test geometry.

The loop and every scattering module read the shared config singleton, so the geometry is
set on it through pytest's monkeypatch (restored after each test). The statistics objects
are no-op stubs: these tests read what they need from the Flight and from the particle.
"""

from math import pi
from types import SimpleNamespace

import numpy as np

from freepaths.config import cf
import freepaths.move
from freepaths.options import ParticleType, SimulationMode
from freepaths.particle import Particle
from freepaths.flight import Flight
import freepaths.run_particle as run_particle_module


SPEED = 1e3
TIMESTEP = 1e-12       # 1 nm per step


def configure(monkeypatch, **overrides):
    """A box W x L x T with all walls, cold side at y = L, absorbing hot side at y = 0"""
    settings = dict(
        width=200e-9, length=1000e-9, thickness=100e-9, timestep=TIMESTEP, number_of_timesteps=200000,
        holes=[], pillars=[], interfaces=[],
        include_right_sidewall=True, include_left_sidewall=True,
        include_top_sidewall=False, include_bottom_sidewall=False,
        cold_side_position_top=True, cold_side_position_bottom=False,
        cold_side_position_right=False, cold_side_position_left=False,
        hot_side_position_top=False, hot_side_position_bottom=True,
        hot_side_position_right=False, hot_side_position_left=False,
        rethermalization_on_hot_sides=False, include_internal_scattering=False,
        low_memory_usage=True, output_path_animation=False, output_scattering_map=False,
        is_two_dimensional_material=False, depletion_width=0.0, max_number_of_scattering_events=None,
        side_wall_roughness=0.0, top_roughness=0.0, hole_roughness=0.0, pillar_roughness=0.0,
        interface_roughness=0.0, is_carrier_electron=True, temp=300.0)
    settings.update(overrides)
    for name, value in settings.items():
        monkeypatch.setattr(cf, name, value, raising=False)
    return cf


class LoopParticle(Particle):
    """A particle for the real loop. An electron always scatters diffusely; a phonon on walls
    of zero roughness always specularly. No internal scattering. Records every position."""

    def __init__(self, position, direction, particle_type=ParticleType.ELECTRON, k=2e9, branch_number=0, f=1e12):
        super().__init__()
        self.type = particle_type
        self.x, self.y, self.z = position
        self.theta = float(np.arctan2(direction[0], direction[1]))
        self.phi = float(np.arcsin(np.clip(direction[2], -1, 1)))
        self.speed = SPEED
        self.k = k
        self.f = f
        self.branch_number = branch_number
        self.time_of_internal_scattering = np.inf
        self.positions = []

    @property
    def wavelength(self):
        return 2 * pi / self.k

    def speed_in(self, material):
        return SPEED

    def assign_frequency(self, material):
        pass

    def assign_internal_scattering_time(self, material):
        self.time_of_internal_scattering = np.inf

    def move(self):
        self.x, self.y, self.z = freepaths.move.move(self, cf.timestep)
        self.positions.append((self.x, self.y, self.z))


def run(particle, monkeypatch, count_diffuse=False):
    """Run the particle through the real loop; returns its Flight (and the number of
    diffuse events the scattering dispatcher reported, if asked)"""
    diffuse_events = [0]
    if count_diffuse:
        original = run_particle_module.surface_scattering

        def counting(pt, scattering_types, places):
            original(pt, scattering_types, places)
            if scattering_types.is_diffuse:
                diffuse_events[0] += 1
        monkeypatch.setattr(run_particle_module, "surface_scattering", counting)

    stub = SimpleNamespace(save_scattering_events=lambda *args: None, record_time_in_segment=lambda *args: None)
    flight = Flight(particle)
    run_particle_module.run_particle(particle, flight, stub, stub, stub, None, None, None, SimulationMode.PHONON_TRACING)
    return (flight, diffuse_events[0]) if count_diffuse else flight
