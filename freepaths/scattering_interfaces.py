"""Module that provides scattering functions on through the interfaces made of other materials"""

from math import pi, cos, sin, acos
from random import random

from freepaths.scattering_types import Scattering

# Single definitions, shared with the free-surface primitives (scatterers.py star-imports
# both modules, so a second definition here would shadow the original):
from freepaths.scattering_primitives import specularity, no_new_scattering, lambert_angles


def vertical_surface_left_scattering_2T(pt, roughness, cf, is_diffuse=False):
    """Transmission through a thin layer (two interfaces) towards the left"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(pt.phi)*sin(abs(pt.theta)))
    p = specularity(a, roughness, pt)

    # Specular transmission: the phonon is refracted on the way into the layer and refracted
    # back on the way out, so it leaves in its original direction:
    if random() < p and (not is_diffuse):
        return Scattering.SPECULAR

    # Diffuse scattering:
    else:

        # Lambert cosine distribution, redrawn if the new direction would immediately leave the domain
        # (near an edge with the ceiling or floor). The event is diffuse whatever the draw:
        for attempt in range(10):
            pt.theta, pt.phi = lambert_angles(-pi/2)
            if no_new_scattering(pt, cf):
                break
        return Scattering.DIFFUSE


def vertical_surface_right_scattering_2T(pt, roughness, cf, is_diffuse=False):
    """Transmission through a thin layer (two interfaces) towards the right"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(pt.phi)*sin(abs(pt.theta)))
    p = specularity(a, roughness, pt)

    # Specular transmission: the phonon is refracted on the way into the layer and refracted
    # back on the way out, so it leaves in its original direction:
    if random() < p and (not is_diffuse):
        return Scattering.SPECULAR

    # Diffuse scattering:
    else:

        # Lambert cosine distribution, redrawn if the new direction would immediately leave the domain
        # (near an edge with the ceiling or floor). The event is diffuse whatever the draw:
        for attempt in range(10):
            pt.theta, pt.phi = lambert_angles(+pi/2)
            if no_new_scattering(pt, cf):
                break
        return Scattering.DIFFUSE
