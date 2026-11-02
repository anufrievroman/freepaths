"""Module that provides scattering functions on basic objects like various walls"""

from math import pi, cos, sin, tan, exp, atan, atan2, asin, acos, sqrt
from random import random
from numpy import sign

from freepaths.options import ParticleType
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


def lambert_angles(normal_theta, normal_phi=0.0):
    """
    Draw a direction from the Lambert cosine law about the surface normal
    given by (normal_theta, normal_phi), in the direction convention of move.py:
    d = (sin(theta)*cos(phi), cos(theta)*cos(phi), sin(phi)).
    The cosine to the normal is mu = sqrt(random()) (density 2*mu), and the azimuth
    around the normal is uniform. Returns the (theta, phi) of the new direction.
    """
    mu = sqrt(random())
    s = sqrt(1 - mu**2)
    psi = 2*pi*random()
    sin_t, cos_t = sin(normal_theta), cos(normal_theta)
    sin_p, cos_p = sin(normal_phi), cos(normal_phi)

    # Normal and two unit tangents (along d/dtheta and d/dphi of the normal):
    d_x = mu*sin_t*cos_p + s*cos(psi)*cos_t - s*sin(psi)*sin_t*sin_p
    d_y = mu*cos_t*cos_p - s*cos(psi)*sin_t - s*sin(psi)*cos_t*sin_p
    d_z = mu*sin_p + s*sin(psi)*cos_p
    return atan2(d_x, d_y), asin(max(-1.0, min(1.0, d_z)))


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
        pt.theta, pt.phi = lambert_angles(-pi/2)

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
        pt.theta, pt.phi = lambert_angles(+pi/2)

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
    pt.theta, pt.phi = lambert_angles(pi)
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
    pt.theta, pt.phi = lambert_angles(0)
    return Scattering.DIFFUSE


def inclined_surfaces_down_scattering(pt, beta, x, x0, roughness):
    """Scattering from a inclined surfaces pointing down like V"""

    # Calculate angle to the surface normal and specular scattering probability:
    normal_theta = pi - sign(x - x0)*(pi/2 - beta)
    a = acos(abs(cos(pt.phi)*cos(pt.theta - normal_theta)))
    p = specularity(a, roughness, pt)

    # Specular scattering:
    if random() < p:
        pt.theta = - pt.theta + sign(x - x0)*2*beta
        return Scattering.SPECULAR

    # Diffuse scattering:
    pt.theta, pt.phi = lambert_angles(normal_theta)
    return Scattering.DIFFUSE


def inclined_surfaces_up_scattering(pt, beta, x, x0, roughness):
    """Scattering from a inclined surfaces pointing up like ^"""

    # Calculate angle to the surface normal and specular scattering probability:
    normal_theta = sign(x - x0)*(pi/2 - beta)
    a = acos(abs(cos(pt.phi)*cos(pt.theta - normal_theta)))
    p = specularity(a, roughness, pt)

    # Specular scattering:
    if random() < p:
        pt.theta = - pt.theta - sign(x - x0)*2*beta
        return Scattering.SPECULAR

    # Diffuse scattering:
    pt.theta, pt.phi = lambert_angles(normal_theta)
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
    pt.theta, pt.phi = lambert_angles(0, -sign(pt.phi)*pi/2)
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
        pt.theta, pt.phi = lambert_angles(tangent_theta - pi*(y < y0))

        # Accept the angles only if they do not immediately cause new scattering:
        if no_new_scattering(pt, cf):
            return Scattering.DIFFUSE


def circle_inner_scattering(pt, tangent_theta, y, y0, roughness):
    """Scattering from the inner surface of the circle"""

    # Calculate angle to the surface normal (radial, at angle tangent_theta up to pi) and specular scattering probability:
    a = acos(abs(cos(pt.phi)*cos(pt.theta - tangent_theta)))
    p = specularity(a, roughness, pt)

    # Specular scattering:
    if random() < p:
        pt.theta = - pt.theta - pi + 2*tangent_theta
        return Scattering.SPECULAR

    # Diffuse scattering:
    pt.theta, pt.phi = lambert_angles(tangent_theta + pi*(y >= y0))
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
