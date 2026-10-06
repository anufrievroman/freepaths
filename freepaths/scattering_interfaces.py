"""Module that provides scattering functions on through the interfaces made of other materials"""

from math import pi, cos, sin, tan, exp, atan, asin, acos
from random import random
from numpy import sign

from freepaths.options import ParticleType
from freepaths.move import move
from freepaths.scattering_types import Scattering
from freepaths.materials import Si, SiGe


def specularity(angle, roughness, particle):
    """Calculate probability of specular scattering with Soffer's equation"""
    if particle.type is ParticleType.ELECTRON: # No Soffer's equation for electrons, scattering is always diffusive
        return 0
    return exp(-16 * pi**2 * roughness**2 * ((cos(angle))**2) / particle.wavelength**2)


def no_new_scattering(pt, cf):
    """
    Check if new angles do not immediately lead to a new top/bottom or sidewall scattering event.
    Such additional scatering event may cause troubles because at this stage we already checked the domain boundaries.
    Thus, this check is necessary to prevent particles leaving the structure boundaries.
    """
    x, y, z = move(pt, cf.timestep)
    return abs(z) < cf.thickness / 2 and abs(x) < cf.width / 2 and cf.length > y > 0


def vertical_surface_left_scattering_2T(pt, roughness, cf, vg_i_to_vg_j, is_diffuse=False):
    """Scattering from a vertical surface to the left"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(pt.phi)*sin(abs(pt.theta)))
    p = specularity(a, roughness, pt)

    theta_i = pi / 2 - pt.theta # because in the paper is the projection angle of the x axis

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

        # Lambert cosine distribution:
        pt.theta = -pi/2 + asin(2*random() - 1)
        pt.phi = asin((asin(2*random() - 1))/(pi/2))

        # Accept the angles if they do not cause new scattering:
        if no_new_scattering(pt, cf):
            return Scattering.DIFFUSE


def vertical_surface_right_scattering_2T(pt, roughness, cf, vg_i_to_vg_j, is_diffuse=False):
    """Scattering from a vertical surface to the left"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(pt.phi)*sin(abs(pt.theta)))
    p = specularity(a, roughness, pt)

    theta_i = pi / 2 - pt.theta # because in the paper is the projection angle of the x axis

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

        # Lambert cosine distribution:
        pt.theta = +pi/2 + asin(2*random() - 1)
        pt.phi = asin((asin(2*random() - 1))/(pi/2))

        # Accept the angles if they do not cause new scattering:
        if no_new_scattering(pt, cf):
            return Scattering.DIFFUSE
