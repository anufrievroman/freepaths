"""Module that provides scattering functions on basic objects like various walls"""

from math import pi, cos, sin, tan, exp, atan, asin, acos
from random import random
from numpy import sign

from freepaths.particle_types import ParticleType
from freepaths.move import move
from freepaths.scattering_types import Scattering


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


def random_scattering(pt):
    """Scattering in a random direction"""
    pt.theta = -pi + random()*2*pi
    pt.phi = asin(2*random() - 1)
    return Scattering.DIFFUSE


def vertical_surface_left_scattering(pt, roughness, cf, is_diffuse=False):
    """Scattering from a vertical surface to the left"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(pt.phi)*sin(abs(pt.theta)))
    p = specularity(a, roughness, pt)

    # Specular scattering:
    if random() < p and (not is_diffuse):
        pt.theta = - pt.theta
        return Scattering.SPECULAR

    # Diffuse scattering
    else:

        # Lambert cosine distribution:
        pt.theta = -pi/2 + asin(2*random() - 1)
        pt.phi = asin((asin(2*random() - 1))/(pi/2))

        # Accept the angles if they do not cause new scattering:
        if no_new_scattering(pt, cf):
            return Scattering.DIFFUSE


def vertical_surface_right_scattering(pt, roughness, cf, is_diffuse=False):
    """Scattering from a vertical surface to the left"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(pt.phi)*sin(abs(pt.theta)))
    p = specularity(a, roughness, pt)

    # Specular scattering:
    if random() < p and (not is_diffuse):
        pt.theta = - pt.theta
        return Scattering.SPECULAR

    # Diffuse scattering:
    else:

        # Lambert cosine distribution:
        pt.theta = +pi/2 + asin(2*random() - 1)
        pt.phi = asin((asin(2*random() - 1))/(pi/2))

        # Accept the angles if they do not cause new scattering:
        if no_new_scattering(pt, cf):
            return Scattering.DIFFUSE


def horizontal_surface_down_scattering(pt, roughness, is_diffuse=False):
    """Scattering from a horizontal surface down"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(pt.phi)*cos(pt.theta))
    p = specularity(a, roughness, pt)

    # Specular scattering:
    if random() < p and (not is_diffuse):
        pt.theta = sign(pt.theta)*pi - pt.theta
        return Scattering.SPECULAR

    # Diffuse scattering:
    # Lambert cosine distribution:
    rand_sign = sign((2*random() - 1))
    pt.theta = rand_sign*pi/2 + rand_sign*acos(random())
    pt.phi = asin((asin(2*random() - 1))/(pi/2))
    return Scattering.DIFFUSE


def horizontal_surface_up_scattering(pt, roughness, is_diffuse=False):
    """Scattering from a horizontal surface up"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(pt.phi)*cos(pt.theta))
    p = specularity(a, roughness, pt)

    # Specular scattering:
    if random() < p and not is_diffuse:
        pt.theta = sign(pt.theta)*pi - pt.theta
        return Scattering.SPECULAR

    # Diffuse scattering:
    # Lambert cosine distribution:
    pt.theta = asin(2*random() - 1)
    pt.phi = asin((asin(2*random() - 1))/(pi/2))
    return Scattering.DIFFUSE


def inclined_surfaces_down_scattering(pt, beta, x, x0, roughness):
    """Scattering from a inclined surfaces pointing down like V"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(pt.phi)*cos(pi/2 - abs(pt.theta) - beta))
    p = specularity(a, roughness, pt)

    # Specular scattering:
    if random() < p:
        pt.theta = - pt.theta + sign(x - x0)*2*beta
        return Scattering.SPECULAR

    # Diffuse scattering:
    rand_sign = sign((2*random() - 1))
    pt.theta = rand_sign*pi - rand_sign*asin(random()) - sign(x - x0)*(pi/2 - beta)
    pt.phi = asin((asin(2*random() - 1))/(pi/2))
    return Scattering.DIFFUSE


def inclined_surfaces_up_scattering(pt, beta, x, x0, roughness):
    """Scattering from a inclined surfaces pointing up like ^"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(pt.phi)*cos(pi/2 - abs(pt.theta) + beta))
    p = specularity(a, roughness, pt)

    # Specular scattering:
    if random() < p:
        pt.theta = - pt.theta - sign(x - x0)*2*beta
        return Scattering.SPECULAR

    # Diffuse scattering:
    pt.theta = asin(2*random() - 1) + sign(x - x0)*(pi/2 - beta)
    pt.phi = asin((asin(2*random() - 1))/(pi/2))
    return Scattering.DIFFUSE


def in_plane_surface_scattering(pt, roughness):
    """Scattering from the ceiling or floor surfaces downwards"""

    # Calculate angle to the surface and specular scattering probability:
    a = pi/2 - abs(pt.phi)
    p = specularity(a, roughness, pt)

    # Specular scattering:
    if random() < p:
        pt.phi = - pt.phi
        return Scattering.SPECULAR

    # Diffuse scattering:
    # Random distribution:
    # theta = -pi + random()*2*pi
    # phi = -asin(random())

    # Lambert cosine distribution:
    pt.theta = - pi + random()*2*pi
    pt.phi = - sign(pt.phi) *  random()*pi/2
    return Scattering.DIFFUSE


def circle_outer_scattering(pt, tangent_theta, y, y0, roughness, cf):
    """Scattering from the outer surface of the circle"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(pt.phi)*cos(pt.theta - tangent_theta))
    p = specularity(a, roughness, pt)

    # Specular scattering:
    if random() < p:
        pt.theta = - pt.theta - pi + 2*tangent_theta
        return Scattering.SPECULAR

    # Diffuse scattering:
    for attempt in range(10):

        # Random distribution:
        # theta = tangent_theta - (-pi/2 + pi*random()) - pi*(y < y0)
        # phi = asin(2*random() - 1)

        # Lambert cosine distribution:
        pt.theta = tangent_theta - (asin(2*random() - 1)) - pi*(y < y0)
        pt.phi = asin((asin(2*random() - 1))/(pi/2))

        # Accept the angles only if they do not immediately cause new scattering:
        if no_new_scattering(pt, cf):
            return Scattering.DIFFUSE


def circle_inner_scattering(pt, tangent_theta, y, y0, roughness):
    """Scattering from the inner surface of the circle"""

    a = atan(tan((pi/2 - pt.theta) + tangent_theta) * cos(pt.phi))
    p = specularity(a, roughness, pt)

    # Specular scattering:
    if random() < p:
        pt.theta = - pt.theta - pi + 2*tangent_theta
        return Scattering.SPECULAR

    # Diffuse scattering:
    pt.theta = tangent_theta - asin(2*random()-1) + pi*(y >= y0)
    pt.phi = asin((asin(2*random() - 1))/(pi/2))
    return Scattering.DIFFUSE


def circle_inclined_inner_scattering(pt, tangent_theta, y, y0, roughness):
    """
    Scattering from the inner inclined surface of the circle.
    This is used for pillars with inclide walls.
    THIS IS NOT TESTED AND NOT WORKING.
    """

    a = atan(tan((pi/2 - pt.theta) + tangent_theta) * cos(pt.phi)) # - (pi / 2 - pillar.wall_angle)))
    p = specularity(a, roughness, pt)

    # Specular scattering:
    if random() < p:

        # If particle moves from the center of the pillar to the wall:
        if distance_from_pillar_center > distance_from_pillar_center_original:

            # If theta does not reflect back:
            if pt.phi < pi/2 - 2 * pillar.wall_angle:
                pt.phi = pt.phi # - (pi / 2 - pillar.wall_angle)

            # Regular reflection:
            else:
                pt.theta = - pt.theta - pi + 2*tangent_theta
                pt.phi = pt.phi # - (pi / 2 - pillar.wall_angle)

        # If particle strikes the wall as it goes towards the center:
        else:
            pt.phi = -sign(pt.phi) * pt.phi #- 2 * pillar.wall_angle
        return Scattering.SPECULAR

    # Diffuse scattering:
    pt.theta = tangent_theta - asin(2*random()-1) + pi*(y >= y0)
    pt.phi = asin((asin(2*random() - 1))/(pi/2)) # - (pi / 2 - pillar.wall_angle)
    return Scattering.DIFFUSE
