"""Module that move a phonon in one timestep using cache of previous moves"""

from numpy import cos, sin
from functools import lru_cache


@lru_cache(maxsize=32)
def step(theta, phi, speed, timestep):
    """Calculate and cache one step of phonon motion"""
    cos_phi = abs(cos(phi))
    d_x = sin(theta) * cos_phi * speed * timestep
    d_y = cos(theta) * cos_phi * speed * timestep
    d_z = sin(phi) * speed * timestep
    return d_x, d_y, d_z


def move(phonon, timestep):
    """Move a phonon in one timestep and return new coordinates"""
    d_x, d_y, d_z = step(phonon.theta, phonon.phi, phonon.speed, timestep)
    new_x = phonon.x + d_x
    new_y = phonon.y + d_y
    new_z = phonon.z + d_z
    return new_x, new_y, new_z
