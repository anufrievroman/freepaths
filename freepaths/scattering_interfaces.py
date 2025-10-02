"""Module that provides scattering functions on through the interfaces made of other materials"""

from math import pi, cos, sin, tan, exp, atan, asin, acos
from random import random
from numpy import sign

from freepaths.particle_types import ParticleType
from freepaths.move import move
from freepaths.scattering_types import Scattering
from freepaths.materials import Si, Ge


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


def vertical_surface_left_scattering_1T(pt, roughness, cf=None, *, mat_in=None, mat_out=None, is_diffuse=False):
    """Scattering from a vertical surface to the left"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(pt.phi)*sin(abs(pt.theta)))
    p = specularity(a, roughness, pt)

    if (mat_in is None) or (mat_out is None):
        mat_0 = next((m for m in cf.materials if getattr(m, "name", "") == "Si"), None) or cf.materials[0]
        mat_j = next((m for m in cf.materials if getattr(m, "name", "") == "Ge"), None) or cf.materials[-1]
    else:
        mat_0, mat_j = mat_in, mat_out

    # Get the original speed of the phonon:
    original_speed = pt.speed
    pt.assign_speed(mat_0);  vg_i = pt.speed
    pt.assign_speed(mat_j); vg_j = pt.speed
    pt.speed = original_speed

    theta_i = pi / 2 - pt.theta # because in the paper is the projection angle of the x axis

    # Specular scattering:
    if random() < p and (not is_diffuse):
        sin_t1 = (vg_i / vg_j) * sin(theta_i)  # theta_i to calcul snell equation with the good angle
        if abs(sin_t1) <= 1: #abs because Snell's law is not defined for sin_t [-1, 1]
            theta_t = -abs(asin(sin_t1))  #  variation of the angle because 1 transmission
        else:
            theta_t = + pt.theta  # total reflexin if impossible angle
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


def vertical_surface_right_scattering_1T(pt, roughness, cf, *, mat_in=None, mat_out=None, is_diffuse=False):
    """Scattering from a vertical surface to the left"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(pt.phi)*sin(abs(pt.theta)))
    p = specularity(a, roughness, pt)

    if (mat_in is None) or (mat_out is None):
        mat_0 = next((m for m in cf.materials if getattr(m, "name", "") == "Si"), None) or cf.materials[0]
        mat_j = next((m for m in cf.materials if getattr(m, "name", "") == "Ge"), None) or cf.materials[-1]
    else:
        mat_0, mat_j = mat_in, mat_out

    # Get the original speed of the phonon:
    original_speed = pt.speed
    pt.assign_speed(mat_0);  vg_i = pt.speed
    pt.assign_speed(mat_j); vg_j = pt.speed
    pt.speed = original_speed

    theta_i = pi / 2 - pt.theta # because in the paper is the projection angle of the x axis

    # Specular scattering:
    if random() < p and (not is_diffuse):
        sin_t1 = (vg_i / vg_j) * sin(theta_i)  # snell law
        if abs(sin_t1) <= 1: #abs because Snell's law is not defined for sin_t [-1, 1]
            theta_t = abs(asin(sin_t1))
        else:
            theta_t = - pt.theta  # Total reflexion if angle impossible
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


def horizontal_surface_up_scattering_1T(pt, roughness, cf, *, mat_in=None, mat_out=None, is_diffuse=False):
    """Interface horizontale: transmission 'vers le haut' (cos(theta)>0) avec SMMM 1T."""

    a = acos(cos(pt.phi) * cos(pt.theta))
    p = specularity(a, roughness, pt)

    if (mat_in is None) or (mat_out is None):
        mat_0 = next((m for m in cf.materials if getattr(m, "name", "") == "Si"), None) or cf.materials[0]
        mat_j = next((m for m in cf.materials if getattr(m, "name", "") == "Ge"), None) or cf.materials[-1]
    else:
        mat_0, mat_j = mat_in, mat_out

    # Get the original speed of the phonon:
    original_speed = pt.speed
    pt.assign_speed(mat_0);  vg_i = pt.speed
    pt.assign_speed(mat_j); vg_j = pt.speed
    pt.speed = original_speed


    theta_i = abs(pt.theta)  # normale = y

    # Specular :
    if (not is_diffuse) and (random() < p):
        sin_t = (vg_i / vg_j) * sin(theta_i)
        if abs(sin_t) <= 1:
            theta_t = abs(asin(sin_t))          # transmitted angle
            pt.theta = sign(pt.theta) * theta_t
        else:
            # Total reflexion horizontal
            pt.theta = sign(pt.theta)*pi - pt.theta
        return Scattering.SPECULAR

    # Diffusive
    pt.theta = asin(2*random() - 1)  # |theta|<=pi/2 => cos(theta)>0
    pt.phi   = asin((asin(2*random() - 1)) / (pi/2))
    return Scattering.DIFFUSE


def horizontal_surface_down_scattering_1T(pt, roughness, cf, *, mat_in=None, mat_out=None, is_diffuse=False):
    """Interface horizontale: transmission 'vers le bas' (cos(theta)<0) avec SMMM 1T."""
    a = acos(cos(pt.phi) * cos(pt.theta))
    p = specularity(a, roughness, pt)

    if (mat_in is None) or (mat_out is None):
        mat_0 = next((m for m in cf.materials if getattr(m, "name", "") == "Si"), None) or cf.materials[0]
        mat_j = next((m for m in cf.materials if getattr(m, "name", "") == "Ge"), None) or cf.materials[-1]
    else:
        mat_0, mat_j = mat_in, mat_out

    # Get the original speed of the phonon:
    original_speed = pt.speed
    pt.assign_speed(mat_0);  vg_i = pt.speed
    pt.assign_speed(mat_j); vg_j = pt.speed
    pt.speed = original_speed
    theta_i = abs(pt.theta)  # normale = y

    print(f"[PRIM][DOWN_1T] in={getattr(mat_0,'name','?')} out={getattr(mat_j,'name','?')} "
          f"vg_i={vg_i:.2e} vg_j={vg_j:.2e}", flush=True)

    # Spéculaire:
    if (not is_diffuse) and (random() < p):
        sin_t = (vg_i / vg_j) * sin(theta_i)
        if abs(sin_t) <= 1:
            theta_t = abs(asin(sin_t))              # (0..pi/2)
            pt.theta = sign(pt.theta)*pi - theta_t  # cos(theta)<0 (to -y)
        else:
            # Total reflexion
            pt.theta = sign(pt.theta)*pi - pt.theta
        return Scattering.SPECULAR

    # Diffusive
    s = sign(2*random() - 1)
    pt.theta = s*pi/2 + s*acos(random())  # cos(theta)<0
    pt.phi   = asin((asin(2*random() - 1)) / (pi/2))
    return Scattering.DIFFUSE


