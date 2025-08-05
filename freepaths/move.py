"""Module that move a particle in one timestep using cache of previous moves"""

from numpy import cos, sin
from functools import lru_cache


@lru_cache(maxsize=32)
def step(theta, phi, speed, timestep):
    """Calculate and cache one step of particle motion"""
    cos_phi = abs(cos(phi))
    d_x = sin(theta) * cos_phi * speed * timestep
    d_y = cos(theta) * cos_phi * speed * timestep
    d_z = sin(phi) * speed * timestep
    return d_x, d_y, d_z


def move(particle, timestep):
    """Move a particle in one timestep and return new coordinates"""
    d_x, d_y, d_z = step(particle.theta, particle.phi, particle.speed, timestep)
    new_x = particle.x + d_x
    new_y = particle.y + d_y
    new_z = particle.z + d_z
    return new_x, new_y, new_z
