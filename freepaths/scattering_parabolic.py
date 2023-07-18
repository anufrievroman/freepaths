"""
Module that provides scattering functions on parabolic boundaries.
The code in this module is contributed by Dhanishta Singh in 2022.
See examples in Singh et al. APL 122, 092203 (2023).
"""

from math import pi, cos, sin, tan, exp, sqrt, atan, asin, acos
from random import random
from numpy import sign

from freepaths.config import cf
from freepaths.move import move
from freepaths.scattering_types import Scattering
from freepaths.scattering_primitives import *
from freepaths.scatterers import *


def top_parabola_scattering(ph, parabola, roughness, scattering_types, x, y, z):
    """Scattering on top parabolic boundary"""

    # If phonon is beyond the parabola:
    y_cept = -(cf.width/2)**2 / (4*parabola.focus) + parabola.tip
    if y > y_cept and (x**2 + 4*parabola.focus*(y - parabola.tip)) >= 0:

        # Calculate angle to the surface and specular scattering probability:
        normal_theta =  pi * (x < 0) - atan(2*parabola.focus/x)
        dot_product = cos(ph.phi) * sin(ph.theta - normal_theta)
        angle = acos(dot_product)
        p = specularity(angle, roughness, ph.wavelength)

        # Specular scattering:
        if random() < p:
            if abs(ph.theta) > pi/2:
                ph.theta = ph.theta - 2*normal_theta
            else :
                ph.theta = 2*normal_theta - ph.theta
            scattering_types.walls = Scattering.SPECULAR

        # Diffuse scattering:
        else :
            scattering_types.walls = Scattering.DIFFUSE
            for attempt in range(10):

                # Lambert distribution
                ph.theta = normal_theta + asin(2*random() - 1) - pi/2
                ph.phi = asin((asin(2*random() - 1))/(pi/2))

                # Accept the angles only if they do not immediately cause new scattering:
                if no_new_scattering(ph):
                    break


def bottom_parabola_scattering(ph, parabola, roughness, scattering_types, x, y, z):
    """Scattering on bottom parabolic boundary"""

    # If phonon is below the parabola:
    y_cept = (cf.width/2)**2 / (4*parabola.focus + parabola.tip)
    if y < y_cept and (x**2 - 4*parabola.focus*(y - parabola.tip)) >= 0:

        # Calculate angle to the surface and specular scattering probability:
        normal_theta =  pi * (x < 0) - atan(-2*parabola.focus/x)
        dot_product = cos(ph.phi) * sin(ph.theta - normal_theta)
        angle = acos(dot_product)
        p = specularity(angle, roughness, ph.wavelength)

        # Specular scattering:
        if random() < p:
            if abs(ph.theta) > pi/2:
                ph.theta = ph.theta - 2*normal_theta
            else :
                ph.theta = 2*normal_theta - ph.theta
            scattering_types.walls = Scattering.SPECULAR

        # Diffuse scattering:
        else :
            scattering_types.walls = Scattering.DIFFUSE
            for attempt in range(10):

                # Lambertian distribution
                ph.theta = normal_theta + asin(2*random() - 1) - pi/2
                ph.phi = asin((asin(2*random() - 1))/(pi/2))

                # Accept the angles only if they do not immediately cause new scattering:
                if no_new_scattering(ph):
                    break
