"""Module that provides scattering functions on through the interfaces made of other materials"""

from math import pi, cos, sin, acos
from random import random

from freepaths.scattering_types import Scattering

# Single definitions, shared with the free-surface primitives. They used to be duplicated
# here verbatim, and since scatterers.py star-imports both modules this copy shadowed the
# original, so any change to the Soffer equation had to be made twice to take effect.
from freepaths.scattering_primitives import specularity, no_new_scattering, lambert_angles


def vertical_surface_left_scattering_2T(pt, roughness, cf, vg_i_to_vg_j, is_diffuse=False):
    """Scattering from a vertical surface to the left"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(pt.phi)*sin(abs(pt.theta)))
    p = specularity(a, roughness, pt)

    # Angle of incidence from the interface normal (x) in 3D. The in-plane pi/2 - theta alone ignored phi,
    # treating out-of-plane phonons as closer to normal incidence than they are:
    theta_i = acos(abs(cos(pt.phi) * sin(pt.theta)))

    # Specular scattering:
    if random() < p and (not is_diffuse):

        sin_t1 = (vg_i_to_vg_j) * sin(theta_i)  # theta_i to calculate snell equation with the good angle
        if abs(sin_t1) <= 1: # abs because Snell's law is not defined for sin_t [-1, 1] 7

            sin_theta_t2 = (1 / vg_i_to_vg_j) * sin_t1
            if abs(sin_theta_t2) <= 1:
                theta_t = + pt.theta
            else:
                theta_t = - pt.theta
        else:
            theta_t = - pt.theta  # totale reflexion if impossible angle

        pt.theta = theta_t
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


def vertical_surface_right_scattering_2T(pt, roughness, cf, vg_i_to_vg_j, is_diffuse=False):
    """Scattering from a vertical surface to the left"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(pt.phi)*sin(abs(pt.theta)))
    p = specularity(a, roughness, pt)

    # Angle of incidence from the interface normal (x) in 3D. The in-plane pi/2 - theta alone ignored phi,
    # treating out-of-plane phonons as closer to normal incidence than they are:
    theta_i = acos(abs(cos(pt.phi) * sin(pt.theta)))

    # Specular scattering:
    if random() < p and (not is_diffuse):
        sin_t1 = (vg_i_to_vg_j) * sin(theta_i)  # snell law
        if abs(sin_t1) <= 1: #abs because Snell's law is not defined for sin_t [-1, 1]
            sin_theta_t2 = (1 / vg_i_to_vg_j) * sin_t1
            if abs(sin_theta_t2) <= 1:
                theta_t = + pt.theta  # no variation of the angle
            else:
                theta_t = - pt.theta
        else:
            theta_t = - pt.theta  # total reflexion if angle impossible
        pt.theta = theta_t
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
