"""Module that provides scattering functions on basic objects like various walls"""

from math import pi, cos, sin, tan, exp, sqrt, atan, asin, acos
from random import random
from numpy import sign

from freepaths.move import move
from freepaths.scattering_types import Scattering


def specularity(angle, roughness, wavelength):
    """Calculate probability of specular scattering with Soffer's equation"""
    return exp(-16 * pi**2 * roughness**2 * ((cos(angle))**2) / wavelength**2)


def no_new_scattering(ph, cf):
    """
    Check if new angles do not immediately lead to a new top/bottom or sidewall scattering event.
    Such additional scatering event may cause troubles because at this stage we already checked the domain boundaries.
    Thus, this check is necessary to prevent phonons leaving the structure boundaries.
    """
    x, y, z = move(ph, cf.timestep)
    return abs(z) < cf.thickness / 2 and abs(x) < cf.width / 2 and cf.length > y > 0


def random_scattering(ph):
    """Scattering in a random direction"""
    ph.theta = -pi + random()*2*pi
    ph.phi = asin(2*random() - 1)
    return Scattering.DIFFUSE


def vertical_surface_left_scattering(ph, roughness, cf, is_diffuse=False):
    """Scattering from a vertical surface to the left"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(ph.phi)*sin(abs(ph.theta)))
    p = specularity(a, roughness, ph.wavelength)

    # Specular scattering:
    if random() < p and (not is_diffuse):
        ph.theta = - ph.theta
        return Scattering.SPECULAR

    # Diffuse scattering:
    for attempt in range(10):

        # Lambert cosine distribution:
        ph.theta = -pi/2 + asin(2*random() - 1)
        ph.phi = asin((asin(2*random() - 1))/(pi/2))

        # Accept the angles if they do not cause new scattering:
        if no_new_scattering(ph, cf):
            return Scattering.DIFFUSE


def vertical_surface_right_scattering(ph, roughness, cf, is_diffuse=False):
    """Scattering from a vertical surface to the left"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(ph.phi)*sin(abs(ph.theta))) # Angle to the surface
    p = specularity(a, roughness, ph.wavelength)

    # Specular scattering:
    if random() < p and (not is_diffuse):
        ph.theta = - ph.theta
        return Scattering.SPECULAR

    # Diffuse scattering:
    for attempt in range(10):

        # Lambert cosine distribution:
        ph.theta = +pi/2 + asin(2*random() - 1)
        ph.phi = asin((asin(2*random() - 1))/(pi/2))

        # Accept the angles if they do not cause new scattering:
        if no_new_scattering(ph, cf):
            return Scattering.DIFFUSE


def horizontal_surface_down_scattering(ph, roughness, is_diffuse=False):
    """Scattering from a horizontal surface down"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(ph.phi)*cos(ph.theta))
    p = specularity(a, roughness, ph.wavelength)

    # Specular scattering:
    if random() < p and (not is_diffuse):
        ph.theta = sign(ph.theta)*pi - ph.theta
        return Scattering.SPECULAR

    # Diffuse scattering:
    # Lambert cosine distribution:
    rand_sign = sign((2*random() - 1))
    ph.theta = rand_sign*pi/2 + rand_sign*acos(random())
    ph.phi = asin((asin(2*random() - 1))/(pi/2))
    return Scattering.DIFFUSE


def horizontal_surface_up_scattering(ph, roughness, is_diffuse=False):
    """Scattering from a horizontal surface up"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(ph.phi)*cos(ph.theta))
    p = specularity(a, roughness, ph.wavelength)

    # Specular scattering:
    if random() < p and not is_diffuse:
        ph.theta = sign(ph.theta)*pi - ph.theta
        return Scattering.SPECULAR

    # Diffuse scattering:
    # Lambert cosine distribution:
    ph.theta = asin(2*random() - 1)
    ph.phi = asin((asin(2*random() - 1))/(pi/2))
    return Scattering.DIFFUSE


def inclined_surfaces_down_scattering(ph, beta, x, x0, roughness):
    """Scattering from a inclined surfaces pointing down like V"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(ph.phi)*cos(pi/2 - abs(ph.theta) - beta))
    p = specularity(a, roughness, ph.wavelength)

    # Specular scattering:
    if random() < p:
        ph.theta = - ph.theta + sign(x - x0)*2*beta
        return Scattering.SPECULAR

    # Diffuse scattering:
    rand_sign = sign((2*random() - 1))
    ph.theta = rand_sign*pi - rand_sign*asin(random()) - sign(x - x0)*(pi/2 - beta)
    ph.phi = asin((asin(2*random() - 1))/(pi/2))
    return Scattering.DIFFUSE


def inclined_surfaces_up_scattering(ph, beta, x, x0, roughness):
    """Scattering from a inclined surfaces pointing up like ^"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(ph.phi)*cos(pi/2 - abs(ph.theta) + beta))
    # a = acos(cos(ph.phi)*cos(abs(ph.theta) - pi/2 + beta))
    p = specularity(a, roughness, ph.wavelength)

    # Specular scattering:
    if random() < p:
        ph.theta = - ph.theta - sign(x - x0)*2*beta
        return Scattering.SPECULAR

    # Diffuse scattering:
    ph.theta = asin(2*random() - 1) + sign(x - x0)*(pi/2 - beta)
    ph.phi = asin((asin(2*random() - 1))/(pi/2))
    return Scattering.DIFFUSE


def in_plane_surface_scattering(ph, roughness):
    """Scattering from the ceiling or floor surfaces downwards"""

    # Calculate angle to the surface and specular scattering probability:
    a = pi/2 - abs(ph.phi)
    p = specularity(a, roughness, ph.wavelength)

    # Specular scattering:
    if random() < p:
        ph.phi = - ph.phi
        return Scattering.SPECULAR

    # Diffuse scattering:
    # Random distribution:
    # theta = -pi + random()*2*pi
    # phi = -asin(random())

    # Lambert cosine distribution:
    ph.theta = - pi + random()*2*pi
    ph.phi = - sign(ph.phi) *  random()*pi/2
    return Scattering.DIFFUSE


def circle_outer_scattering(ph, tangent_theta, y, y0, roughness, cf):
    """Scattering from the outer surface of the circle"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(ph.phi)*cos(ph.theta - tangent_theta))
    p = specularity(a, roughness, ph.wavelength)

    # Specular scattering:
    if random() < p:
        ph.theta = - ph.theta - pi + 2*tangent_theta
        return Scattering.SPECULAR

    # Diffuse scattering:
    for attempt in range(10):

        # Random distribution:
        # theta = tangent_theta - (-pi/2 + pi*random()) - pi*(y < y0)
        # phi = asin(2*random() - 1)

        # Lambert cosine distribution:
        ph.theta = tangent_theta - (asin(2*random() - 1)) - pi*(y < y0)
        ph.phi = asin((asin(2*random() - 1))/(pi/2))

        # Accept the angles only if they do not immediately cause new scattering:
        if no_new_scattering(ph, cf):
            return Scattering.DIFFUSE


def circle_inner_scattering(ph, tangent_theta, y, y0, roughness):
    """Scattering from the inner surface of the circle"""

    a = atan(tan((pi/2 - ph.theta) + tangent_theta) * cos(ph.phi))
    p = specularity(a, roughness, ph.wavelength)

    # Specular scattering:
    if random() < p:
        ph.theta = - ph.theta - pi + 2*tangent_theta
        return Scattering.SPECULAR

    # Diffuse scattering:
    ph.theta = tangent_theta - asin(2*random()-1) + pi*(y >= y0)
    ph.phi = asin((asin(2*random() - 1))/(pi/2))
    return Scattering.DIFFUSE


def circle_inclined_inner_scattering(ph, tangent_theta, y, y0, roughness):
    """
    Scattering from the inner inclined surface of the circle.
    This is used for pillars with inclide walls.
    THIS IS NOT TESTED AND NOT WORKING.
    """

    a = atan(tan((pi/2 - ph.theta) + tangent_theta) * cos(ph.phi)) # - (pi / 2 - pillar.wall_angle)))
    p = specularity(a, roughness, ph.wavelength)

    # Specular scattering:
    if random() < p:

        # If phonon moves from the center of the pillar to the wall:
        if distance_from_pillar_center > distance_from_pillar_center_original:

            # If theta does not reflect back:
            if ph.phi < pi/2 - 2 * pillar.wall_angle:
                ph.phi = ph.phi # - (pi / 2 - pillar.wall_angle)

            # Regular reflection:
            else:
                ph.theta = - ph.theta - pi + 2*tangent_theta
                ph.phi = ph.phi # - (pi / 2 - pillar.wall_angle)

        # If phonon strikes the wall as it goes towards the center:
        else:
            ph.phi = -sign(ph.phi) * ph.phi #- 2 * pillar.wall_angle
        return Scattering.SPECULAR

    # Diffuse scattering:
    ph.theta = tangent_theta - asin(2*random()-1) + pi*(y >= y0)
    ph.phi = asin((asin(2*random() - 1))/(pi/2)) # - (pi / 2 - pillar.wall_angle)
    return Scattering.DIFFUSE
