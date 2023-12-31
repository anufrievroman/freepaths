"""
Module that provides scattering functions on parabolic semicircular holes.
The code in this module is contributed by Felix Barbier in 2023.
"""

from math import pi, cos, sin, tan, sqrt, atan, asin, acos
from random import random
from numpy import sign

from freepaths.scattering_types import Scattering
from freepaths.scattering_primitives import *


def scattering_on_semicircular_holes(ph, x0, y0, R, scattering_types, x, y, z, cf):
    """Check if a phonon strikes a circular hole and calculate the new direction"""

    # If phonon is inside the circle with radius R:
    if (x - x0) ** 2 + (y - y0) ** 2 <= R**2 and x >= x0:
        r1 = abs(y0 - y + (x - x0) * cos(ph.theta) / abs(sin(ph.theta)))
        if R > r1:
            # Calculate angle to the surface and specular scattering probability:
            a = acos(
                cos(ph.phi) * sin(abs(ph.theta))
            )  # Angle to the normal to the surface
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                scattering_types.holes = Scattering.SPECULAR
                ph.theta = -ph.theta

            # Diffuse scattering:
            else:
                scattering_types.holes = Scattering.DIFFUSE
                attempt = 0
                while attempt < 10:
                    attempt += 1
                    # Lambert cosine distribution:
                    ph.theta = -sign(sin(ph.theta)) * pi / 2 + asin(2 * random() - 1)
                    ph.phi = asin((asin(2 * random() - 1)) / (pi / 2))

                    # Accept the angles only if they do not lead to new scattering:
                    if no_new_scattering(ph, cf):
                        break
        else:
            # Calculate angle to the surface and specular scattering probability:
            if y == y0:
                y += 1e-9  # Prevent division by zero
            tangent_theta = atan((x - x0) / (y - y0))
            a = acos(cos(ph.phi) * cos(ph.theta + sign(y - y0) * tangent_theta))
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = -ph.theta - pi + 2 * tangent_theta
                scattering_types.holes = Scattering.SPECULAR

            # Diffuse scattering:
            else:
                scattering_types.holes = Scattering.DIFFUSE
                attempt = 0
                while attempt < 10:
                    attempt += 1

                    # Random distribution:
                    # theta = tangent_theta - (-pi/2 + pi*random()) - pi*(y < y0)
                    # phi = asin(2*random() - 1)

                    # Lambert cosine distribution:
                    ph.theta = tangent_theta - (asin(2 * random() - 1)) - pi * (y < y0)
                    ph.phi = asin((asin(2 * random() - 1)) / (pi / 2))

                    # Accept the angles only if they do not immediately cause new scattering:
                    if no_new_scattering(ph, cf):
                        break


def scattering_on_arccircular_v_holes(
    ph, x0, y0, R, Rinner, alphap, scattering_types, x, y, z, cf
):
    """Check if a phonon strikes a circular hole and calculate the new direction"""
    if x == x0:  # to prevent division by 0
        x = x + 1e-12
    theta0 = atan((y - y0) / (x - x0))
    if y == y0:  # to prevent division by 0
        y = y + 1e-12
    tangent_theta = atan(
        (x - x0) / (y - y0)
    )  # use for scatering sepcular of circular boundary

    # to know if it is in or not
    if (
        (Rinner**2 <= (x - x0) ** 2 + (y - y0) ** 2 <= R**2)
        and x >= x0
        and -alphap / 2 <= theta0 <= alphap / 2
    ):
        xp = ph.x - x0
        yp = ph.y - y0
        thetapre = atan((yp) / (xp))  # angle of previous point

        # Parameter
        continu = 0  # to know if already scatered
        beta = pi / 2 - alphap / 2  # angle for triangle
        pillar_wall_angle = pi / 2

        if (
            yp > 0
            and (R**2 <= xp**2 + yp**2 and -alphap / 2 <= thetapre <= alphap / 2)
            == False
            and (
                Rinner**2 >= xp**2 + yp**2 and -alphap / 2 < thetapre < alphap / 2
            )
            == False
        ):
            continu = 1
            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi) * cos(ph.theta - (pi / 2 - beta)))
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = -ph.theta + 2 * beta
                scattering_types.holes = Scattering.SPECULAR

            else:
                ph.theta = asin(2 * random() - 1) - 1 * (pi / 2 - beta)
                ph.phi = asin((asin(2 * random() - 1)) / (pi / 2))
                scattering_types.holes = Scattering.DIFFUSE

        if (
            yp < 0
            and (R**2 <= xp**2 + yp**2 and -alphap / 2 <= thetapre <= alphap / 2)
            == False
            and (
                Rinner**2 >= xp**2 + yp**2 and -alphap / 2 < thetapre < alphap / 2
            )
            == False
        ):
            continu = 1

            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi) * cos(ph.theta + (pi / 2 - beta)))
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = -ph.theta - 2 * beta
                scattering_types.holes = Scattering.SPECULAR
            else:
                # rand_sign = sign((2*random() - 1))
                ph.theta = pi - asin(random()) + 1 * (pi / 2 - beta)
                ph.phi = asin((asin(2 * random() - 1)) / (pi / 2))
                scattering_types.holes = Scattering.DIFFUSE

        if R**2 <= xp**2 + yp**2 and continu == 0:
            continu = 1
            if y == y0:
                y += 1e-9  # Prevent division by zero
            tangent_theta = atan((x - x0) / (y - y0))
            a = acos(cos(ph.phi) * cos(ph.theta + sign(y - y0) * tangent_theta))
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = -ph.theta - pi + 2 * tangent_theta
                scattering_types.holes = Scattering.SPECULAR
            else:
                scattering_types.holes = Scattering.DIFFUSE
                attempt = 0
                while attempt < 10:
                    attempt += 1

                    # Random distribution:
                    # theta = tangent_theta - (-pi/2 + pi*random()) - pi*(y < y0)
                    # phi = asin(2*random() - 1)
                    # #Lambert cosine distribution:
                    ph.theta = tangent_theta - (asin(2 * random() - 1)) - pi * (y < y0)
                    ph.phi = asin((asin(2 * random() - 1)) / (pi / 2))

                    # Accept the angles only if they do not immediately cause new scattering:
                    if no_new_scattering(ph, cf):
                        break

        if continu == 0:
            tangent_theta = atan((x - x0) / (y - y0))
            a = atan(
                tan((pi / 2 - ph.theta) + tangent_theta)
                * cos(ph.phi - (pi / 2 - pillar_wall_angle))
            )
            p = specularity(a, cf.pillar_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                # If phonon moves from the center of the pillar to the wall:
                if sqrt((abs(x) - abs(x0)) ** 2 + (abs(y) - abs(y0)) ** 2) >= sqrt(
                    (abs(ph.x) - abs(x0)) ** 2 + (abs(ph.y) - abs(y0)) ** 2
                ):
                    # If theta does not reflect back:
                    # if ph.phi < pi/2 - 2 * pillar_wall_angle:
                    # ph.phi = ph.phi - (pi / 2 - pillar_wall_angle)

                    # Regular reflection:
                    # else:
                    ph.theta = -ph.theta - pi + 2 * tangent_theta
                    ph.phi = ph.phi - (pi / 2 - pillar_wall_angle)

                # If phonon strikes the wall as it goes towards the center:
                else:
                    ph.theta = -ph.theta - pi + 2 * tangent_theta
                    ph.phi = -sign(ph.phi) * ph.phi - 2 * pillar_wall_angle
                scattering_types.holes = Scattering.SPECULAR

            else:
                # attempt = 0
                # while attempt < 10:
                # attempt += 1

                ph.theta = tangent_theta - asin(2 * random() - 1) + pi * (y >= y0)
                ph.phi = asin((asin(2 * random() - 1)) / (pi / 2)) - (
                    pi / 2 - cf.pillar_wall_angle
                )
                scattering_types.pillars = Scattering.DIFFUSE
                # if no_new_scattering(ph, cf):
                # break


def scattering_on_arccircular_v_demi_down_holes(
    ph, x0, y0, R, Rinner, alphap, alphap2, scattering_types, x, y, z, cf
):
    """Check if a phonon strikes a circular hole and calculate the new direction"""
    if x == x0:  # to prevent division by 0
        x = x + 1e-12
    theta0 = atan((y - y0) / (x - x0))
    if y == y0:  # to prevent division by 0
        y = y + 1e-12
    tangent_theta = atan(
        (x - x0) / (y - y0)
    )  # use for scatering sepcular of circular boundary

    # to know if it is in or not
    if (
        (Rinner**2 <= (x - x0) ** 2 + (y - y0) ** 2 <= R**2)
        and x >= x0
        and -alphap / 2 <= theta0 <= -alphap2 / 2
    ):
        xp = ph.x - x0
        yp = ph.y - y0
        thetapre = atan((yp) / (xp))  # angle of previous point

        # Parameter
        continu = 0  # to know if already scatered
        beta = pi / 2 - alphap / 2  # angle for triangle
        beta2 = pi / 2 - alphap2 / 2
        pillar_wall_angle = pi / 2
        y_mid = -xp * tan(alphap / 2) + 10e-9  # -alphap2/2)

        if (
            yp > y_mid
            and (
                R**2 <= xp**2 + yp**2 and -alphap / 2 <= thetapre <= -alphap2 / 2
            )
            == False
            and (
                Rinner**2 >= xp**2 + yp**2
                and -alphap / 2 < thetapre < -alphap2 / 2
            )
            == False
        ):
            continu = 1
            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi) * cos(ph.theta + (pi / 2 - beta2)))
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = -ph.theta - 2 * beta2
                scattering_types.holes = Scattering.SPECULAR

            else:
                ph.theta = asin(2 * random() - 1) + 1 * (pi / 2 - beta2)
                ph.phi = asin((asin(2 * random() - 1)) / (pi / 2))
                scattering_types.holes = Scattering.DIFFUSE

        if (
            yp < y_mid
            and (
                R**2 <= xp**2 + yp**2 and -alphap / 2 <= thetapre <= -alphap2 / 2
            )
            == False
            and (
                Rinner**2 >= xp**2 + yp**2
                and -alphap / 2 < thetapre < -alphap2 / 2
            )
            == False
        ):
            continu = 1

            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi) * cos(ph.theta + (pi / 2 - beta)))
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = -ph.theta - 2 * beta
                scattering_types.holes = Scattering.SPECULAR
            else:
                # rand_sign = sign((2*random() - 1))
                ph.theta = pi - asin(random()) + 1 * (pi / 2 - beta)
                ph.phi = asin((asin(2 * random() - 1)) / (pi / 2))
                scattering_types.holes = Scattering.DIFFUSE

        if R**2 <= xp**2 + yp**2 and continu == 0:
            # for sol in S:
            # Ima_sol = sol- conj(sol)
            # if Ima_sol==0:
            # if -1<=sol<=1:
            # if abs(acos(sol))<=alphap/2:

            continu = 1
            if y == y0:
                y += 1e-9  # Prevent division by zero
            tangent_theta = atan((x - x0) / (y - y0))
            a = acos(cos(ph.phi) * cos(ph.theta + sign(y - y0) * tangent_theta))
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = -ph.theta - pi + 2 * tangent_theta
                scattering_types.holes = Scattering.SPECULAR
            else:
                scattering_types.holes = Scattering.DIFFUSE
                attempt = 0
                while attempt < 10:
                    attempt += 1

                    # Random distribution:
                    # theta = tangent_theta - (-pi/2 + pi*random()) - pi*(y < y0)
                    # phi = asin(2*random() - 1)
                    # #Lambert cosine distribution:
                    ph.theta = tangent_theta - (asin(2 * random() - 1)) - pi * (y < y0)
                    ph.phi = asin((asin(2 * random() - 1)) / (pi / 2))

                    # Accept the angles only if they do not immediately cause new scattering:
                    if no_new_scattering(ph, cf):
                        break

        if continu == 0:
            tangent_theta = atan((x - x0) / (y - y0))
            a = atan(
                tan((pi / 2 - ph.theta) + tangent_theta)
                * cos(ph.phi - (pi / 2 - pillar_wall_angle))
            )
            p = specularity(a, cf.pillar_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                # If phonon moves from the center of the pillar to the wall:
                if sqrt((abs(x) - abs(x0)) ** 2 + (abs(y) - abs(y0)) ** 2) >= sqrt(
                    (abs(ph.x) - abs(x0)) ** 2 + (abs(ph.y) - abs(y0)) ** 2
                ):
                    # If theta does not reflect back:
                    # if ph.phi < pi/2 - 2 * pillar_wall_angle:
                    # ph.phi = ph.phi - (pi / 2 - pillar_wall_angle)

                    # Regular reflection:
                    # else:
                    ph.theta = -ph.theta - pi + 2 * tangent_theta
                    ph.phi = ph.phi - (pi / 2 - pillar_wall_angle)

                # If phonon strikes the wall as it goes towards the center:
                else:
                    ph.theta = -ph.theta - pi + 2 * tangent_theta
                    ph.phi = -sign(ph.phi) * ph.phi - 2 * pillar_wall_angle
                scattering_types.holes = Scattering.SPECULAR

            else:
                # attempt = 0
                # while attempt < 10:
                # attempt += 1

                ph.theta = tangent_theta - asin(2 * random() - 1) + pi * (y >= y0)
                ph.phi = asin((asin(2 * random() - 1)) / (pi / 2)) - (
                    pi / 2 - cf.pillar_wall_angle
                )
                scattering_types.pillars = Scattering.DIFFUSE
                # if no_new_scattering(ph, cf):


def scattering_on_arccircular_v_demi_up_holes(
    ph, x0, y0, R, Rinner, alphap, alphap2, scattering_types, x, y, z, cf
):
    """Check if a phonon strikes a circular hole and calculate the new direction"""
    if x == x0:  # to prevent division by 0
        x = x + 1e-12
    theta0 = atan((y - y0) / (x - x0))
    if y == y0:  # to prevent division by 0
        y = y + 1e-12
    tangent_theta = atan(
        (x - x0) / (y - y0)
    )  # use for scatering sepcular of circular boundary

    # to know if it is in or not
    if (
        (Rinner**2 <= (x - x0) ** 2 + (y - y0) ** 2 <= R**2)
        and x >= x0
        and alphap2 / 2 <= theta0 <= alphap / 2
    ):
        xp = ph.x - x0
        yp = ph.y - y0
        thetapre = atan((yp) / (xp))  # angle of previous point

        # Parameter
        continu = 0  # to know if already scatered
        beta = pi / 2 - alphap / 2  # angle for triangle
        beta2 = pi / 2 - alphap2 / 2
        pillar_wall_angle = pi / 2
        y_mid = xp * tan(alphap / 2) - 10e-9  # -alphap2/2)

        if (
            yp > y_mid
            and (R**2 <= xp**2 + yp**2 and alphap2 / 2 <= thetapre <= alphap / 2)
            == False
            and (
                Rinner**2 >= xp**2 + yp**2 and alphap2 / 2 < thetapre < alphap / 2
            )
            == False
        ):
            continu = 1
            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi) * cos(ph.theta - (pi / 2 - beta)))
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = -ph.theta + 2 * beta
                scattering_types.holes = Scattering.SPECULAR

            else:
                ph.theta = asin(2 * random() - 1) - 1 * (pi / 2 - beta)
                ph.phi = asin((asin(2 * random() - 1)) / (pi / 2))
                scattering_types.holes = Scattering.DIFFUSE

        if (
            yp < y_mid
            and (R**2 <= xp**2 + yp**2 and alphap2 / 2 <= thetapre <= alphap / 2)
            == False
            and (
                Rinner**2 >= xp**2 + yp**2 and alphap2 / 2 < thetapre < alphap / 2
            )
            == False
        ):
            continu = 1

            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi) * cos(ph.theta - (pi / 2 - beta2)))
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = -ph.theta + 2 * beta2
                scattering_types.holes = Scattering.SPECULAR
            else:
                # rand_sign = sign((2*random() - 1))
                ph.theta = pi - asin(random()) - 1 * (pi / 2 - beta2)
                ph.phi = asin((asin(2 * random() - 1)) / (pi / 2))
                scattering_types.holes = Scattering.DIFFUSE

        if R**2 <= xp**2 + yp**2 and continu == 0:
            # for sol in S:
            # Ima_sol = sol- conj(sol)
            # if Ima_sol==0:
            # if -1<=sol<=1:
            # if abs(acos(sol))<=alphap/2:

            continu = 1
            if y == y0:
                y += 1e-9  # Prevent division by zero
            tangent_theta = atan((x - x0) / (y - y0))
            a = acos(cos(ph.phi) * cos(ph.theta + sign(y - y0) * tangent_theta))
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = -ph.theta - pi + 2 * tangent_theta
                scattering_types.holes = Scattering.SPECULAR
            else:
                scattering_types.holes = Scattering.DIFFUSE
                attempt = 0
                while attempt < 10:
                    attempt += 1

                    # Random distribution:
                    # theta = tangent_theta - (-pi/2 + pi*random()) - pi*(y < y0)
                    # phi = asin(2*random() - 1)
                    # #Lambert cosine distribution:
                    ph.theta = tangent_theta - (asin(2 * random() - 1)) - pi * (y < y0)
                    ph.phi = asin((asin(2 * random() - 1)) / (pi / 2))

                    # Accept the angles only if they do not immediately cause new scattering:
                    if no_new_scattering(ph, cf):
                        break

        if continu == 0:
            tangent_theta = atan((x - x0) / (y - y0))
            a = atan(
                tan((pi / 2 - ph.theta) + tangent_theta)
                * cos(ph.phi - (pi / 2 - pillar_wall_angle))
            )
            p = specularity(a, cf.pillar_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                # If phonon moves from the center of the pillar to the wall:
                if sqrt((abs(x) - abs(x0)) ** 2 + (abs(y) - abs(y0)) ** 2) >= sqrt(
                    (abs(ph.x) - abs(x0)) ** 2 + (abs(ph.y) - abs(y0)) ** 2
                ):
                    # If theta does not reflect back:
                    # if ph.phi < pi/2 - 2 * pillar_wall_angle:
                    # ph.phi = ph.phi - (pi / 2 - pillar_wall_angle)

                    # Regular reflection:
                    # else:
                    ph.theta = -ph.theta - pi + 2 * tangent_theta
                    ph.phi = ph.phi - (pi / 2 - pillar_wall_angle)

                # If phonon strikes the wall as it goes towards the center:
                else:
                    ph.theta = -ph.theta - pi + 2 * tangent_theta
                    ph.phi = -sign(ph.phi) * ph.phi - 2 * pillar_wall_angle
                scattering_types.holes = Scattering.SPECULAR

            else:
                # attempt = 0
                # while attempt < 10:
                # attempt += 1

                ph.theta = tangent_theta - asin(2 * random() - 1) + pi * (y >= y0)
                ph.phi = asin((asin(2 * random() - 1)) / (pi / 2)) - (
                    pi / 2 - cf.pillar_wall_angle
                )
                scattering_types.pillars = Scattering.DIFFUSE
                # if no_new_scattering(ph, cf):
                # break                                 #break


def scattering_on_arccircular_h_holes(
    ph, x0, y0, R, Rinner, alphap, scattering_types, x, y, z, cf
):
    """Check if a phonon strikes a circular hole and calculate the new direction"""
    if x == x0:  # to prevent division by 0
        x = x + 1e-12
    theta0 = atan((x - x0) / (y - y0))
    if y == y0:  # to prevent division by 0
        y = y + 1e-12
    tangent_theta = atan((x - x0) / (y - y0))

    # to know if it is in or not
    if (
        (Rinner**2 <= (x - x0) ** 2 + (y - y0) ** 2 <= R**2)
        and y >= y0
        and -alphap / 2 <= theta0 <= alphap / 2
    ):
        xp = ph.x - x0
        yp = ph.y - y0
        thetapre = atan((xp) / (yp))

        continu = 0
        beta = alphap / 2
        pillar_wall_angle = pi / 2

        if (
            xp > 0
            and (R**2 <= xp**2 + yp**2 and -alphap / 2 <= thetapre <= alphap / 2)
            == False
            and (
                Rinner**2 >= xp**2 + yp**2 and -alphap / 2 < thetapre < alphap / 2
            )
            == False
        ):
            continu = 1
            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi) * cos(ph.theta - (pi / 2 - beta)))
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = -ph.theta + 2 * beta
                scattering_types.holes = Scattering.SPECULAR

            else:
                ph.theta = pi - (asin(2 * random() - 1) + 1 * (pi / 2 - beta))
                ph.phi = asin((asin(2 * random() - 1)) / (pi / 2))
                scattering_types.holes = Scattering.DIFFUSE

        if (
            xp < 0
            and (R**2 <= xp**2 + yp**2 and -alphap / 2 <= thetapre <= alphap / 2)
            == False
            and (
                Rinner**2 >= xp**2 + yp**2 and -alphap / 2 < thetapre < alphap / 2
            )
            == False
        ):
            continu = 1

            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi) * cos(ph.theta + (pi / 2 - beta)))
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = -ph.theta - 2 * beta
                scattering_types.holes = Scattering.SPECULAR
            else:
                # rand_sign = sign((2*random() - 1))
                ph.theta = pi - asin(random()) + 1 * (pi / 2 - beta)
                ph.phi = asin((asin(2 * random() - 1)) / (pi / 2))
                scattering_types.holes = Scattering.DIFFUSE

        if R**2 <= xp**2 + yp**2 and continu == 0:
            # for sol in S:
            # Ima_sol = sol- conj(sol)
            # if Ima_sol==0:
            # if -1<=sol<=1:
            # if abs(acos(sol))<=alphap/2:

            continu = 1
            if y == y0:
                y += 1e-9  # Prevent division by zero
            tangent_theta = atan((x - x0) / (y - y0))
            a = acos(cos(ph.phi) * cos(ph.theta + sign(y - y0) * tangent_theta))
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = -ph.theta - pi + 2 * tangent_theta
                scattering_types.holes = Scattering.SPECULAR
            else:
                scattering_types.holes = Scattering.DIFFUSE
                attempt = 0
                while attempt < 10:
                    attempt += 1

                    # Random distribution:
                    # theta = tangent_theta - (-pi/2 + pi*random()) - pi*(y < y0)
                    # phi = asin(2*random() - 1)
                    # #Lambert cosine distribution:
                    ph.theta = tangent_theta - (asin(2 * random() - 1)) - pi * (y < y0)
                    ph.phi = asin((asin(2 * random() - 1)) / (pi / 2))

                    # Accept the angles only if they do not immediately cause new scattering:
                    if no_new_scattering(ph, cf):
                        break

        if continu == 0:
            tangent_theta = atan((x - x0) / (y - y0))
            a = atan(
                tan((pi / 2 - ph.theta) + tangent_theta)
                * cos(ph.phi - (pi / 2 - pillar_wall_angle))
            )
            p = specularity(a, cf.pillar_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                # If phonon moves from the center of the pillar to the wall:
                if sqrt((abs(x) - abs(x0)) ** 2 + (abs(y) - abs(y0)) ** 2) >= sqrt(
                    (abs(ph.x) - abs(x0)) ** 2 + (abs(ph.y) - abs(y0)) ** 2
                ):
                    # If theta does not reflect back:
                    # if ph.phi < pi/2 - 2 * pillar_wall_angle:
                    # ph.phi = ph.phi - (pi / 2 - pillar_wall_angle)

                    # Regular reflection:
                    # else:
                    ph.theta = -ph.theta - pi + 2 * tangent_theta
                    ph.phi = ph.phi - (pi / 2 - pillar_wall_angle)

                # If phonon strikes the wall as it goes towards the center:
                else:
                    ph.theta = -ph.theta - pi + 2 * tangent_theta
                    ph.phi = -sign(ph.phi) * ph.phi - 2 * pillar_wall_angle
                scattering_types.holes = Scattering.SPECULAR

            else:
                # attempt = 0
                # while attempt < 10:
                # attempt += 1

                ph.theta = tangent_theta - asin(2 * random() - 1) + pi * (y >= y0)
                ph.phi = asin((asin(2 * random() - 1)) / (pi / 2)) - (
                    pi / 2 - cf.pillar_wall_angle
                )
                scattering_types.pillars = Scattering.DIFFUSE
                # if no_new_scattering(ph, cf):


def scattering_on_arccircular_h_reverse_holes(
    ph, x0, y0, R, Rinner, alphap, scattering_types, x, y, z, cf
):
    """Check if a phonon strikes a circular hole and calculate the new direction"""
    if x == x0:  # to prevent division by 0
        x = x + 1e-12
    theta0 = atan((x - x0) / (y - y0))
    if y == y0:  # to prevent division by 0
        y = y + 1e-12
    tangent_theta = atan((x - x0) / (y - y0))

    # to know if it is in or not
    if (
        (Rinner**2 <= (x - x0) ** 2 + (y - y0) ** 2 <= R**2)
        and y <= y0
        and -alphap / 2 <= theta0 <= alphap / 2
    ):
        xp = ph.x - x0
        yp = ph.y - y0
        thetapre = atan((xp) / (yp))

        continu = 0
        beta = alphap / 2
        pillar_wall_angle = pi / 2

        if (
            xp > 0
            and (R**2 <= xp**2 + yp**2 and -alphap / 2 <= thetapre <= alphap / 2)
            == False
            and (
                Rinner**2 >= xp**2 + yp**2 and -alphap / 2 < thetapre < alphap / 2
            )
            == False
        ):
            continu = 1
            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi) * cos(ph.theta + (pi / 2 - beta)))
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = -ph.theta - 2 * beta
                scattering_types.holes = Scattering.SPECULAR

            else:
                ph.theta = asin(2 * random() - 1) + 1 * (pi / 2 - beta)
                ph.phi = asin((asin(2 * random() - 1)) / (pi / 2))
                scattering_types.holes = Scattering.DIFFUSE

        if (
            xp < 0
            and (R**2 <= xp**2 + yp**2 and -alphap / 2 <= thetapre <= alphap / 2)
            == False
            and (
                Rinner**2 >= xp**2 + yp**2 and -alphap / 2 < thetapre < alphap / 2
            )
            == False
        ):
            continu = 1

            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi) * cos(ph.theta - (pi / 2 - beta)))
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = -ph.theta + 2 * beta
                scattering_types.holes = Scattering.SPECULAR
            else:
                # rand_sign = sign((2*random() - 1))

                ph.theta = asin(2 * random() - 1) - 1 * (pi / 2 - beta)
                ph.phi = asin((asin(2 * random() - 1)) / (pi / 2))
                scattering_types.holes = Scattering.DIFFUSE

        if R**2 <= xp**2 + yp**2 and continu == 0:
            # for sol in S:
            # Ima_sol = sol- conj(sol)
            # if Ima_sol==0:
            # if -1<=sol<=1:
            # if abs(acos(sol))<=alphap/2:

            continu = 1
            if y == y0:
                y += 1e-9  # Prevent division by zero
            tangent_theta = atan((x - x0) / (y - y0))
            a = acos(cos(ph.phi) * cos(ph.theta + sign(y - y0) * tangent_theta))
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = -ph.theta - pi + 2 * tangent_theta
                scattering_types.holes = Scattering.SPECULAR
            else:
                scattering_types.holes = Scattering.DIFFUSE
                attempt = 0
                while attempt < 10:
                    attempt += 1

                    # Random distribution:
                    # theta = tangent_theta - (-pi/2 + pi*random()) - pi*(y < y0)
                    # phi = asin(2*random() - 1)
                    # #Lambert cosine distribution:
                    ph.theta = tangent_theta - (asin(2 * random() - 1)) - pi * (y < y0)
                    ph.phi = asin((asin(2 * random() - 1)) / (pi / 2))

                    # Accept the angles only if they do not immediately cause new scattering:
                    if no_new_scattering(ph, cf):
                        break

        if continu == 0:
            tangent_theta = atan((x - x0) / (y - y0))
            a = atan(
                tan((pi / 2 - ph.theta) + tangent_theta)
                * cos(ph.phi - (pi / 2 - pillar_wall_angle))
            )
            p = specularity(a, cf.pillar_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                # If phonon moves from the center of the pillar to the wall:
                if sqrt((abs(x) - abs(x0)) ** 2 + (abs(y) - abs(y0)) ** 2) >= sqrt(
                    (abs(ph.x) - abs(x0)) ** 2 + (abs(ph.y) - abs(y0)) ** 2
                ):
                    # If theta does not reflect back:
                    # if ph.phi < pi/2 - 2 * pillar_wall_angle:
                    # ph.phi = ph.phi - (pi / 2 - pillar_wall_angle)

                    # Regular reflection:
                    # else:
                    ph.theta = -ph.theta - pi + 2 * tangent_theta
                    ph.phi = ph.phi - (pi / 2 - pillar_wall_angle)

                # If phonon strikes the wall as it goes towards the center:
                else:
                    ph.theta = -ph.theta - pi + 2 * tangent_theta
                    ph.phi = -sign(ph.phi) * ph.phi - 2 * pillar_wall_angle
                scattering_types.holes = Scattering.SPECULAR

            else:
                # attempt = 0
                # while attempt < 10:
                # attempt += 1

                ph.theta = tangent_theta - asin(2 * random() - 1) + pi * (y >= y0)
                ph.phi = asin((asin(2 * random() - 1)) / (pi / 2)) - (
                    pi / 2 - cf.pillar_wall_angle
                )
                scattering_types.pillars = Scattering.DIFFUSE
                # if no_new_scattering(ph, cf):
                # break
