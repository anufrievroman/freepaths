from math import cos, sin
from functools import lru_cache

from parameters import timestep

@lru_cache(maxsize=32)
def step(theta, phi, speed):
    '''Calculate and cache one step of phonon motion'''
    cos_phi = abs(cos(phi))
    dx = sin(theta)*cos_phi*speed*timestep
    dy = cos(theta)*cos_phi*speed*timestep
    dz = sin(phi)*speed*timestep
    return dx, dy, dz


def move(x, y, z, theta, phi, speed):
    '''This function moves a phonon in one timestep and returns new coordinates'''
    dx, dy, dz = step(theta, phi, speed)
    x += dx
    y += dy
    z += dz
    return x, y, z
