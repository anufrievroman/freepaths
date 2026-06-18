"""
Module that provides scattering functions on parabolic semicircular holes.
The code in this module is contributed by Felix Barbier in 2023.
"""

from math import pi, cos, sin, tan, sqrt, atan, asin, acos
from random import random
from numpy import sign

from freepaths.scattering_types import Scattering
from freepaths.scattering_primitives import *


def _scatter_on_circular_arc(pt, x, y, x0, y0, scattering_types, cf):
    """Handle scattering on the circular arc boundary of an arc-shaped hole."""
    if y == y0:
        y += 1e-9
    tangent_theta = atan((x - x0) / (y - y0))
    a = acos(cos(pt.phi) * cos(pt.theta + sign(y - y0) * tangent_theta))
    p = specularity(a, cf.hole_roughness, pt)
    if random() < p:
        pt.theta = -pt.theta - pi + 2 * tangent_theta
        scattering_types.holes = Scattering.SPECULAR
    else:
        scattering_types.holes = Scattering.DIFFUSE
        for _ in range(10):
            pt.theta = tangent_theta - asin(2 * random() - 1) - pi * (y < y0)
            pt.phi = asin((asin(2 * random() - 1)) / (pi / 2))
            if no_new_scattering(pt, cf):
                break


def _scatter_on_pillar_wall(pt, x, y, x0, y0, pillar_wall_angle, scattering_types, cf):
    """Handle scattering on the flat pillar wall when the particle is inside the arc region."""
    tangent_theta = atan((x - x0) / (y - y0))
    a = atan(
        tan((pi / 2 - pt.theta) + tangent_theta)
        * cos(pt.phi - (pi / 2 - pillar_wall_angle))
    )
    p = specularity(a, cf.pillar_roughness, pt)
    if random() < p:
        if sqrt((abs(x) - abs(x0)) ** 2 + (abs(y) - abs(y0)) ** 2) >= sqrt(
            (abs(pt.x) - abs(x0)) ** 2 + (abs(pt.y) - abs(y0)) ** 2
        ):
            pt.theta = -pt.theta - pi + 2 * tangent_theta
            pt.phi = pt.phi - (pi / 2 - pillar_wall_angle)
        else:
            pt.theta = -pt.theta - pi + 2 * tangent_theta
            pt.phi = -sign(pt.phi) * pt.phi - 2 * pillar_wall_angle
        scattering_types.holes = Scattering.SPECULAR
    else:
        pt.theta = tangent_theta - asin(2 * random() - 1) + pi * (y >= y0)
        pt.phi = asin((asin(2 * random() - 1)) / (pi / 2)) - (pi / 2 - pillar_wall_angle)
        scattering_types.pillars = Scattering.DIFFUSE


def scattering_on_semicircular_holes(pt, x0, y0, R, scattering_types, x, y, z, cf):
    """Check if a particle strikes a circular hole and calculate the new direction"""

    if (x - x0) ** 2 + (y - y0) ** 2 <= R**2 and x >= x0:
        r1 = abs(y0 - y + (x - x0) * cos(pt.theta) / abs(sin(pt.theta)))
        if R > r1:
            a = acos(cos(pt.phi) * sin(abs(pt.theta)))
            p = specularity(a, cf.hole_roughness, pt)

            # Specular scattering:
            if random() < p:
                scattering_types.holes = Scattering.SPECULAR
                pt.theta = -pt.theta

            # Diffuse scattering:
            else:
                scattering_types.holes = Scattering.DIFFUSE
                for _ in range(10):
                    pt.theta = -sign(sin(pt.theta)) * pi / 2 + asin(2 * random() - 1)
                    pt.phi = asin((asin(2 * random() - 1)) / (pi / 2))
                    if no_new_scattering(pt, cf):
                        break
        else:
            if y == y0:
                y += 1e-9
            tangent_theta = atan((x - x0) / (y - y0))
            a = acos(cos(pt.phi) * cos(pt.theta + sign(y - y0) * tangent_theta))
            p = specularity(a, cf.hole_roughness, pt)

            # Specular scattering:
            if random() < p:
                pt.theta = -pt.theta - pi + 2 * tangent_theta
                scattering_types.holes = Scattering.SPECULAR

            # Diffuse scattering:
            else:
                scattering_types.holes = Scattering.DIFFUSE
                for _ in range(10):
                    pt.theta = tangent_theta - asin(2 * random() - 1) - pi * (y < y0)
                    pt.phi = asin((asin(2 * random() - 1)) / (pi / 2))
                    if no_new_scattering(pt, cf):
                        break


def scattering_on_arccircular_v_holes(
    pt, x0, y0, R, Rinner, alphap, scattering_types, x, y, z, cf
):
    """Check if a particle strikes a vertical arc-shaped hole and calculate the new direction"""
    if x == x0:
        x += 1e-12
    theta0 = atan((y - y0) / (x - x0))
    if y == y0:
        y += 1e-12
    tangent_theta = atan((x - x0) / (y - y0))

    if (
        (Rinner**2 <= (x - x0) ** 2 + (y - y0) ** 2 <= R**2)
        and x >= x0
        and -alphap / 2 <= theta0 <= alphap / 2
    ):
        xp = pt.x - x0
        yp = pt.y - y0
        thetapre = atan(yp / xp)

        continu = 0
        beta = pi / 2 - alphap / 2
        pillar_wall_angle = pi / 2

        if (
            yp > 0
            and not (R**2 <= xp**2 + yp**2 and -alphap / 2 <= thetapre <= alphap / 2)
            and not (Rinner**2 >= xp**2 + yp**2 and -alphap / 2 < thetapre < alphap / 2)
        ):
            continu = 1
            a = acos(cos(pt.phi) * cos(pt.theta - (pi / 2 - beta)))
            p = specularity(a, cf.hole_roughness, pt)
            if random() < p:
                pt.theta = -pt.theta + 2 * beta
                scattering_types.holes = Scattering.SPECULAR
            else:
                pt.theta = asin(2 * random() - 1) - 1 * (pi / 2 - beta)
                pt.phi = asin((asin(2 * random() - 1)) / (pi / 2))
                scattering_types.holes = Scattering.DIFFUSE

        if (
            yp < 0
            and not (R**2 <= xp**2 + yp**2 and -alphap / 2 <= thetapre <= alphap / 2)
            and not (Rinner**2 >= xp**2 + yp**2 and -alphap / 2 < thetapre < alphap / 2)
        ):
            continu = 1
            a = acos(cos(pt.phi) * cos(pt.theta + (pi / 2 - beta)))
            p = specularity(a, cf.hole_roughness, pt)
            if random() < p:
                pt.theta = -pt.theta - 2 * beta
                scattering_types.holes = Scattering.SPECULAR
            else:
                pt.theta = pi - asin(random()) + 1 * (pi / 2 - beta)
                pt.phi = asin((asin(2 * random() - 1)) / (pi / 2))
                scattering_types.holes = Scattering.DIFFUSE

        if R**2 <= xp**2 + yp**2 and continu == 0:
            continu = 1
            _scatter_on_circular_arc(pt, x, y, x0, y0, scattering_types, cf)

        if continu == 0:
            _scatter_on_pillar_wall(pt, x, y, x0, y0, pillar_wall_angle, scattering_types, cf)


def scattering_on_arccircular_v_demi_down_holes(
    pt, x0, y0, R, Rinner, alphap, alphap2, scattering_types, x, y, z, cf
):
    """Check if a particle strikes a downward vertical demi-arc hole and calculate the new direction"""
    if x == x0:
        x += 1e-12
    theta0 = atan((y - y0) / (x - x0))
    if y == y0:
        y += 1e-12
    tangent_theta = atan((x - x0) / (y - y0))

    if (
        (Rinner**2 <= (x - x0) ** 2 + (y - y0) ** 2 <= R**2)
        and x >= x0
        and -alphap / 2 <= theta0 <= -alphap2 / 2
    ):
        xp = pt.x - x0
        yp = pt.y - y0
        thetapre = atan(yp / xp)

        continu = 0
        beta = pi / 2 - alphap / 2
        beta2 = pi / 2 - alphap2 / 2
        pillar_wall_angle = pi / 2
        y_mid = -xp * tan(alphap / 2) + 10e-9

        if (
            yp > y_mid
            and not (R**2 <= xp**2 + yp**2 and -alphap / 2 <= thetapre <= -alphap2 / 2)
            and not (Rinner**2 >= xp**2 + yp**2 and -alphap / 2 < thetapre < -alphap2 / 2)
        ):
            continu = 1
            a = acos(cos(pt.phi) * cos(pt.theta + (pi / 2 - beta2)))
            p = specularity(a, cf.hole_roughness, pt)
            if random() < p:
                pt.theta = -pt.theta - 2 * beta2
                scattering_types.holes = Scattering.SPECULAR
            else:
                pt.theta = asin(2 * random() - 1) + 1 * (pi / 2 - beta2)
                pt.phi = asin((asin(2 * random() - 1)) / (pi / 2))
                scattering_types.holes = Scattering.DIFFUSE

        if (
            yp < y_mid
            and not (R**2 <= xp**2 + yp**2 and -alphap / 2 <= thetapre <= -alphap2 / 2)
            and not (Rinner**2 >= xp**2 + yp**2 and -alphap / 2 < thetapre < -alphap2 / 2)
        ):
            continu = 1
            a = acos(cos(pt.phi) * cos(pt.theta + (pi / 2 - beta)))
            p = specularity(a, cf.hole_roughness, pt)
            if random() < p:
                pt.theta = -pt.theta - 2 * beta
                scattering_types.holes = Scattering.SPECULAR
            else:
                pt.theta = pi - asin(random()) + 1 * (pi / 2 - beta)
                pt.phi = asin((asin(2 * random() - 1)) / (pi / 2))
                scattering_types.holes = Scattering.DIFFUSE

        if R**2 <= xp**2 + yp**2 and continu == 0:
            continu = 1
            _scatter_on_circular_arc(pt, x, y, x0, y0, scattering_types, cf)

        if continu == 0:
            _scatter_on_pillar_wall(pt, x, y, x0, y0, pillar_wall_angle, scattering_types, cf)


def scattering_on_arccircular_v_demi_up_holes(
    pt, x0, y0, R, Rinner, alphap, alphap2, scattering_types, x, y, z, cf
):
    """Check if a particle strikes an upward vertical demi-arc hole and calculate the new direction"""
    if x == x0:
        x += 1e-12
    theta0 = atan((y - y0) / (x - x0))
    if y == y0:
        y += 1e-12
    tangent_theta = atan((x - x0) / (y - y0))

    if (
        (Rinner**2 <= (x - x0) ** 2 + (y - y0) ** 2 <= R**2)
        and x >= x0
        and alphap2 / 2 <= theta0 <= alphap / 2
    ):
        xp = pt.x - x0
        yp = pt.y - y0
        thetapre = atan(yp / xp)

        continu = 0
        beta = pi / 2 - alphap / 2
        beta2 = pi / 2 - alphap2 / 2
        pillar_wall_angle = pi / 2
        y_mid = xp * tan(alphap / 2) - 10e-9

        if (
            yp > y_mid
            and not (R**2 <= xp**2 + yp**2 and alphap2 / 2 <= thetapre <= alphap / 2)
            and not (Rinner**2 >= xp**2 + yp**2 and alphap2 / 2 < thetapre < alphap / 2)
        ):
            continu = 1
            a = acos(cos(pt.phi) * cos(pt.theta - (pi / 2 - beta)))
            p = specularity(a, cf.hole_roughness, pt)
            if random() < p:
                pt.theta = -pt.theta + 2 * beta
                scattering_types.holes = Scattering.SPECULAR
            else:
                pt.theta = asin(2 * random() - 1) - 1 * (pi / 2 - beta)
                pt.phi = asin((asin(2 * random() - 1)) / (pi / 2))
                scattering_types.holes = Scattering.DIFFUSE

        if (
            yp < y_mid
            and not (R**2 <= xp**2 + yp**2 and alphap2 / 2 <= thetapre <= alphap / 2)
            and not (Rinner**2 >= xp**2 + yp**2 and alphap2 / 2 < thetapre < alphap / 2)
        ):
            continu = 1
            a = acos(cos(pt.phi) * cos(pt.theta - (pi / 2 - beta2)))
            p = specularity(a, cf.hole_roughness, pt)
            if random() < p:
                pt.theta = -pt.theta + 2 * beta2
                scattering_types.holes = Scattering.SPECULAR
            else:
                pt.theta = pi - asin(random()) - 1 * (pi / 2 - beta2)
                pt.phi = asin((asin(2 * random() - 1)) / (pi / 2))
                scattering_types.holes = Scattering.DIFFUSE

        if R**2 <= xp**2 + yp**2 and continu == 0:
            continu = 1
            _scatter_on_circular_arc(pt, x, y, x0, y0, scattering_types, cf)

        if continu == 0:
            _scatter_on_pillar_wall(pt, x, y, x0, y0, pillar_wall_angle, scattering_types, cf)


def scattering_on_arccircular_h_holes(
    pt, x0, y0, R, Rinner, alphap, scattering_types, x, y, z, cf
):
    """Check if a particle strikes a horizontal arc-shaped hole and calculate the new direction"""
    if x == x0:
        x += 1e-12
    theta0 = atan((x - x0) / (y - y0))
    if y == y0:
        y += 1e-12
    tangent_theta = atan((x - x0) / (y - y0))

    if (
        (Rinner**2 <= (x - x0) ** 2 + (y - y0) ** 2 <= R**2)
        and y >= y0
        and -alphap / 2 <= theta0 <= alphap / 2
    ):
        xp = pt.x - x0
        yp = pt.y - y0
        thetapre = atan(xp / yp)

        continu = 0
        beta = alphap / 2
        pillar_wall_angle = pi / 2

        if (
            xp > 0
            and not (R**2 <= xp**2 + yp**2 and -alphap / 2 <= thetapre <= alphap / 2)
            and not (Rinner**2 >= xp**2 + yp**2 and -alphap / 2 < thetapre < alphap / 2)
        ):
            continu = 1
            a = acos(cos(pt.phi) * cos(pt.theta - (pi / 2 - beta)))
            p = specularity(a, cf.hole_roughness, pt)
            if random() < p:
                pt.theta = -pt.theta + 2 * beta
                scattering_types.holes = Scattering.SPECULAR
            else:
                pt.theta = pi - (asin(2 * random() - 1) + 1 * (pi / 2 - beta))
                pt.phi = asin((asin(2 * random() - 1)) / (pi / 2))
                scattering_types.holes = Scattering.DIFFUSE

        if (
            xp < 0
            and not (R**2 <= xp**2 + yp**2 and -alphap / 2 <= thetapre <= alphap / 2)
            and not (Rinner**2 >= xp**2 + yp**2 and -alphap / 2 < thetapre < alphap / 2)
        ):
            continu = 1
            a = acos(cos(pt.phi) * cos(pt.theta + (pi / 2 - beta)))
            p = specularity(a, cf.hole_roughness, pt)
            if random() < p:
                pt.theta = -pt.theta - 2 * beta
                scattering_types.holes = Scattering.SPECULAR
            else:
                pt.theta = pi - asin(random()) + 1 * (pi / 2 - beta)
                pt.phi = asin((asin(2 * random() - 1)) / (pi / 2))
                scattering_types.holes = Scattering.DIFFUSE

        if R**2 <= xp**2 + yp**2 and continu == 0:
            continu = 1
            _scatter_on_circular_arc(pt, x, y, x0, y0, scattering_types, cf)

        if continu == 0:
            _scatter_on_pillar_wall(pt, x, y, x0, y0, pillar_wall_angle, scattering_types, cf)


def scattering_on_arccircular_h_reverse_holes(
    pt, x0, y0, R, Rinner, alphap, scattering_types, x, y, z, cf
):
    """Check if a particle strikes a reversed horizontal arc-shaped hole and calculate the new direction"""
    if x == x0:
        x += 1e-12
    theta0 = atan((x - x0) / (y - y0))
    if y == y0:
        y += 1e-12
    tangent_theta = atan((x - x0) / (y - y0))

    if (
        (Rinner**2 <= (x - x0) ** 2 + (y - y0) ** 2 <= R**2)
        and y <= y0
        and -alphap / 2 <= theta0 <= alphap / 2
    ):
        xp = pt.x - x0
        yp = pt.y - y0
        thetapre = atan(xp / yp)

        continu = 0
        beta = alphap / 2
        pillar_wall_angle = pi / 2

        if (
            xp > 0
            and not (R**2 <= xp**2 + yp**2 and -alphap / 2 <= thetapre <= alphap / 2)
            and not (Rinner**2 >= xp**2 + yp**2 and -alphap / 2 < thetapre < alphap / 2)
        ):
            continu = 1
            a = acos(cos(pt.phi) * cos(pt.theta + (pi / 2 - beta)))
            p = specularity(a, cf.hole_roughness, pt)
            if random() < p:
                pt.theta = -pt.theta - 2 * beta
                scattering_types.holes = Scattering.SPECULAR
            else:
                pt.theta = asin(2 * random() - 1) + 1 * (pi / 2 - beta)
                pt.phi = asin((asin(2 * random() - 1)) / (pi / 2))
                scattering_types.holes = Scattering.DIFFUSE

        if (
            xp < 0
            and not (R**2 <= xp**2 + yp**2 and -alphap / 2 <= thetapre <= alphap / 2)
            and not (Rinner**2 >= xp**2 + yp**2 and -alphap / 2 < thetapre < alphap / 2)
        ):
            continu = 1
            a = acos(cos(pt.phi) * cos(pt.theta - (pi / 2 - beta)))
            p = specularity(a, cf.hole_roughness, pt)
            if random() < p:
                pt.theta = -pt.theta + 2 * beta
                scattering_types.holes = Scattering.SPECULAR
            else:
                pt.theta = asin(2 * random() - 1) - 1 * (pi / 2 - beta)
                pt.phi = asin((asin(2 * random() - 1)) / (pi / 2))
                scattering_types.holes = Scattering.DIFFUSE

        if R**2 <= xp**2 + yp**2 and continu == 0:
            continu = 1
            _scatter_on_circular_arc(pt, x, y, x0, y0, scattering_types, cf)

        if continu == 0:
            _scatter_on_pillar_wall(pt, x, y, x0, y0, pillar_wall_angle, scattering_types, cf)
