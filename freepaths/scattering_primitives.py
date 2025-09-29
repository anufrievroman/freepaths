"""Module that provides scattering functions on basic objects like various walls"""

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


def random_scattering(pt):
    """Scattering in a random direction"""
    pt.theta = -pi + random()*2*pi
    pt.phi = asin(2*random() - 1)
    return Scattering.DIFFUSE



def vertical_surface_left_scattering_2T(pt, roughness, cf, is_diffuse=False):
    """Scattering from a vertical surface to the left"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(pt.phi)*sin(abs(pt.theta)))
    p = specularity(a, roughness, pt)
    

    # Get the group velocities of the materials:
    mat_0 = next((mat for mat in cf.materials if mat.name == "Si"), None)      
    mat_j = next((mat for mat in cf.materials if mat.name == "Ge"), None)       
    
    # Get the original speed of the phonon:
    original_speed = pt.speed
    pt.assign_speed(mat_0)
    vg_i = pt.speed

    pt.assign_speed(mat_j)
    vg_j = pt.speed
    pt.speed = original_speed

    theta_i = pi / 2 - pt.theta # because in the paper is the projection angle of the x axis


    # Specular scattering:
    if random() < p and (not is_diffuse):
       

        sin_t1 = (vg_i / vg_j) * sin(theta_i)  # theta_i to calcul snell equation with the good angle
        if abs(sin_t1) <= 1: #abs because Snell's law is not defined for sin_t [-1, 1] 7
            
            sin_theta_t2 = (vg_j / vg_i) * sin_t1
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


def vertical_surface_right_scattering_2T(pt, roughness, cf, is_diffuse=False):
    """Scattering from a vertical surface to the left"""

    # Calculate angle to the surface and specular scattering probability:
    a = acos(cos(pt.phi)*sin(abs(pt.theta)))
    p = specularity(a, roughness, pt)
    


    # Get the group velocities of the materials:
    mat_0 = next((mat for mat in cf.materials if mat.name == "Si"), None)       # Left material (ex: Si)
    mat_j = next((mat for mat in cf.materials if mat.name == "Ge"), None)       # current materail (ex: Ge)
    
    # Get the original speed of the phonon:
    original_speed = pt.speed
    pt.assign_speed(mat_0)
    vg_i = pt.speed

    pt.assign_speed(mat_j)
    vg_j = pt.speed
    pt.speed = original_speed

    theta_i = pi / 2 - pt.theta # because in the paper is the projection angle of the x axis

    # Specular scattering:
    if random() < p and (not is_diffuse):
        sin_t1 = (vg_i / vg_j) * sin(theta_i)  # snell law
        if abs(sin_t1) <= 1: #abs because Snell's law is not defined for sin_t [-1, 1] 
            sin_theta_t2 = (vg_j / vg_i) * sin_t1
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

    # Diffuse scattering
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
