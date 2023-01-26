"""Modules provides scattering processes on various objects"""

from math import pi, cos, sin, tan, exp, sqrt, atan, asin, acos
from random import random
from numpy import sign

from freepaths.config import cf
from freepaths.move import move
from freepaths.scattering_types import Scattering


def specularity(angle, roughness, wavelength):
    """Calculate probability of specular scattering with Soffer's equation"""
    return exp(-16 * pi**2 * roughness**2 * ((cos(angle))**2) / wavelength**2)


def internal_scattering(ph, flight, scattering_types):
    """This function is checking if the time passed since previous diffuse scattering event
    reached the time until an internal scattering event, and if yes, scatters randomly"""
    if flight.time_since_previous_scattering >= ph.time_of_internal_scattering:
        ph.theta = -pi + random()*2*pi
        ph.phi = asin(2*random() - 1)
        scattering_types.internal = Scattering.DIFFUSE


def reinitialization(ph, scattering_types):
    """Re-thermalize phonon if it comes back to the hot side"""
    _, y, _ = move(ph, cf.timestep)

    # If phonon returns to the staring line y = 0, generate it again:
    if y < 0:
        ph.assign_angles()
        scattering_types.hot_side = Scattering.DIFFUSE


def top_parabola_scattering(ph, scattering_types):
    """Scattering on top parabolic boundary"""
    x, y, z = move(ph, cf.timestep)

    # If phonon is beyond the parabola:
    y_cept = -(cf.width/2)**2 / (4*cf.top_parabola_focus) + cf.top_parabola_tip
    if y > y_cept and (x**2 + 4*cf.top_parabola_focus*(y - cf.top_parabola_tip)) >= 0:

        # Calculate angle to the surface and specular scattering probability:
        normal_theta =  pi * (x < 0) - atan(2*cf.top_parabola_focus/x)
        dot_product = cos(ph.phi) * sin(ph.theta - normal_theta)
        angle = acos(dot_product)
        p = specularity(angle, cf.side_wall_roughness, ph.wavelength)

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
            attempt = 0
            while attempt < 10:
                attempt += 1

                # Lambertian distribution
                ph.theta = normal_theta + asin(2*random() - 1) - pi/2
                ph.phi = asin((asin(2*random() - 1))/(pi/2))

                # Accept the angles only if they do not immediately cause new scattering:
                if no_new_scattering(ph):
                    break


def bottom_parabola_scattering(ph, scattering_types):
    """Scattering on bottom parabolic boundary"""
    x, y, z = move(ph, cf.timestep)

    # If phonon is below the parabola:
    y_cept = (cf.width/2)**2 / (4*cf.bottom_parabola_focus + cf.bottom_parabola_tip)
    if y < y_cept and (x**2 - 4*cf.bottom_parabola_focus*(y - cf.bottom_parabola_tip)) >= 0:

        # Calculate angle to the surface and specular scattering probability:
        normal_theta =  pi * (x < 0) - atan(-2*cf.bottom_parabola_focus/x)
        dot_product = cos(ph.phi) * sin(ph.theta - normal_theta)
        angle = acos(dot_product)
        p = specularity(angle, cf.side_wall_roughness, ph.wavelength)

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
            attempt = 0
            while attempt < 10:
                attempt += 1

                # Lambertian distribution
                ph.theta = normal_theta + asin(2*random() - 1) - pi/2
                ph.phi = asin((asin(2*random() - 1))/(pi/2))

                # Accept the angles only if they do not immediately cause new scattering:
                if no_new_scattering(ph):
                    break


def scattering_on_circular_holes(ph, x0, y0, R, scattering_types, x, y, z):
    """Check if a phonon strikes a circular hole and calculate the new direction"""

    # If phonon is inside the circle with radius R:
    if (x - x0)**2 + (y - y0)**2 <= R**2:

        # Calculate angle to the surface and specular scattering probability:
        if y == y0: y += 1e-9 # Prevent division by zero
        tangent_theta = atan((x - x0)/(y - y0))
        a = acos(cos(ph.phi)*cos(ph.theta + sign(y - y0)*tangent_theta))
        p = specularity(a, cf.hole_roughness, ph.wavelength)

        # Specular scattering:
        if random() < p:
            ph.theta = - ph.theta - pi + 2*tangent_theta
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
                ph.theta = tangent_theta - (asin(2*random() - 1)) - pi*(y < y0)
                ph.phi = asin((asin(2*random() - 1))/(pi/2))

                # Accept the angles only if they do not immediately cause new scattering:
                if no_new_scattering(ph):
                    break


def scattering_on_rectangular_holes(ph, x0, y0, Lx, Ly, scattering_types, x, y, z):
    """Check if the phonon strikes a rectangular hole and calculate new direction"""

    # If the phonon is inside the rectangle:
    if (abs(x - x0) <= Lx / 2) and (abs(y - y0) <= Ly / 2):

        # Coordinate y of the intersection with the hole side:
        y1 = (y0 - y) + cos(ph.theta)*(Lx/2 - abs(x0 - x))/abs(sin(ph.theta))

        # Scattering on left and right walls of the hole:
        if abs(y1) <= Ly/2:

            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi)*sin(abs(ph.theta)))  # Angle to the normal to the surface
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                scattering_types.holes = Scattering.SPECULAR
                ph.theta = - ph.theta

            # Diffuse scattering:
            else:
                scattering_types.holes = Scattering.DIFFUSE
                attempt = 0
                while attempt < 10:
                    attempt += 1

                    # Lambert cosine distribution:
                    ph.theta = - sign(sin(ph.theta))*pi/2 + asin(2*random()-1)
                    ph.phi = asin((asin(2*random() - 1))/(pi/2))

                    # Accept the angles only if they do not lead to new scattering:
                    if no_new_scattering(ph):
                        break

        # Scattering on top and bottom walls of the hole:
        else:
            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi)*cos(ph.theta))
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                scattering_types.holes = Scattering.SPECULAR
                ph.theta = sign(ph.theta)*pi - ph.theta

            # Diffuse scattering:
            else:
                scattering_types.holes = Scattering.DIFFUSE
                attempt = 0
                while attempt < 10:
                    attempt += 1

                    # Scattering on the top surface of the hole:
                    if abs(ph.theta) > pi / 2:
                        # Lambert cosine distribution:
                        ph.theta = asin(2*random() - 1)
                        ph.phi = asin((asin(2*random() - 1))/(pi/2))

                    # Scattering on the bottom surface of the hole:
                    else:
                        # Lambert cosine distribution:
                        rand_sign = sign((2*random() - 1))
                        ph.theta = rand_sign*pi/2 + rand_sign*acos(random())
                        ph.phi = asin((asin(2*random() - 1))/(pi/2))

                    # Accept the angles only if they do not immediately cause new scattering:
                    if no_new_scattering(ph):
                        break


def scattering_on_circular_pillars(ph, x0, y0, R_base, scattering_types, x, y, z):
    """Check if a phonon strikes a circular pillar and calculate new direction"""

    # Cone radius at a given z coordinate:
    R = R_base - (z - cf.thickness / 2) / tan(cf.pillar_wall_angle)

    # If phonon crosses the pillar boundary. Third condition is to exclude all other pillars:
    if ((x - x0)**2 + (y - y0)**2 >= R**2) and (z > cf.thickness / 2) \
            and ((x - x0) ** 2 + (y - y0) ** 2 < (R + 2 * ph.speed * cf.timestep) ** 2):

        # Calculate angle to the surface and specular scattering probability:
        tangent_theta = atan((x - x0)/(y - y0))
        a = atan(tan((pi/2 - ph.theta) + tangent_theta) * cos(ph.phi - (pi / 2 - cf.pillar_wall_angle)))
        p = specularity(a, cf.pillar_roughness, ph.wavelength)

        # Specular scattering:
        if random() < p:

            # If phonon moves from the center of the pillar to the wall:
            if sqrt((abs(x) - abs(x0))**2 + (abs(y) - abs(y0))**2) >= sqrt((abs(ph.x) - abs(x0))**2 + (abs(ph.y) - abs(y0))**2) :

                # If theta does not reflect back:
                if ph.phi < pi/2 - 2 * cf.pillar_wall_angle:
                    ph.phi = ph.phi - (pi / 2 - cf.pillar_wall_angle)

                # Regular reflection:
                else:
                    ph.theta = - ph.theta - pi + 2*tangent_theta
                    ph.phi = ph.phi - (pi / 2 - cf.pillar_wall_angle)

            # If phonon strikes the wall as it goes towards the center:
            else:
                ph.phi = -sign(ph.phi) * ph.phi - 2 * cf.pillar_wall_angle
            scattering_types.pillars = Scattering.SPECULAR

        # Diffuse scattering:
        else:
            ph.theta = tangent_theta - asin(2*random()-1) + pi*(y >= y0)
            ph.phi = asin((asin(2*random() - 1))/(pi/2)) - (pi / 2 - cf.pillar_wall_angle)
            scattering_types.pillars = Scattering.DIFFUSE


def scattering_on_triangle_down_holes(ph, x0, y0, Lx, Ly, scattering_types, x, y, z):
    """Check if the phonon strikes a reverse triangular hole and calculate new direction after the scattering"""

    # Angle of the triangle:
    beta = atan(0.5*Lx/Ly)

    # If phonon is inside the triangle:
    if (Ly/2 - (y - y0) <= (Lx/2 - abs(x - x0))/tan(beta)) and (abs(y - y0) < Ly/2):

        # Scattering on the bottom wall of the triangle:
        if (y + cf.timestep * ph.speed > y0 + Ly / 2) and (abs(ph.theta) > pi / 2):

            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi)*cos(ph.theta))
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = sign(ph.theta)*pi - ph.theta
                scattering_types.holes = Scattering.SPECULAR

            # Diffuse scattering:
            else:
                # Lambert cosine distribution:
                ph.theta = asin(2*random() - 1)
                ph.phi = asin((asin(2*random() - 1))/(pi/2))
                scattering_types.holes = Scattering.DIFFUSE

        # Scattering on the sidewalls of the triangle:
        else:
            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi)*cos(ph.theta - sign(x - x0)*(pi/2 - beta)))
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = - ph.theta + sign(x - x0)*2*beta
                scattering_types.holes = Scattering.SPECULAR

            # Diffuse scattering:
            else:
                rand_sign = sign((2*random() - 1))
                ph.theta = rand_sign*pi - rand_sign*asin(random()) - sign(x-x0)*(pi/2 - beta)
                ph.phi = asin((asin(2*random() - 1))/(pi/2))
                scattering_types.holes= Scattering.DIFFUSE


def scattering_on_triangle_up_holes(ph, x0, y0, Lx, Ly, scattering_types, x, y, z):
    """Check if the phonon strikes a reverse triangular hole and calculate new direction"""

    # Angle of the triangle:
    beta = atan(0.5*Lx/Ly)

    # If phonon is inside the triangle:
    if (Ly/2 + (y - y0) <= (Lx/2 - abs(x - x0))/tan(beta)) and (abs(y - y0) < Ly/2):
        # x1=Ly/2/tan(theta) - abs(y0-y)/tan(theta) + abs(x0-x)

        # Scattering on the bottom wall of the triangle:
        if ((y - cf.timestep * ph.speed) < (y0 - Ly / 2)) and (abs(ph.theta) < pi / 2):

            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi)*cos(ph.theta))
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = sign(ph.theta)*pi - ph.theta
                scattering_types.holes = Scattering.SPECULAR

            # Diffuse scattering:
            else:
                # Lambert cosine distribution:
                rand_sign = sign((2*random() - 1))
                ph.theta = rand_sign*pi/2 + rand_sign*acos(random())
                ph.phi = asin((asin(2*random() - 1))/(pi/2))
                scattering_types.holes = Scattering.DIFFUSE

        # Scattering on the sidewalls of the triangle:
        else:

            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi)*cos(ph.theta + sign(x - x0)*(pi/2 - beta)))
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = - ph.theta - sign(x - x0)*2*beta
                scattering_types.holes = Scattering.SPECULAR

            # Diffuse scattering:
            else:
                # Lambert cosine distribution:
                ph.theta = asin(2*random() - 1) + sign(x - x0)*(pi/2 - beta)
                ph.phi = asin((asin(2*random() - 1))/(pi/2))
                scattering_types.holes = Scattering.DIFFUSE


def no_new_scattering(ph):
    """Check if new angles do not immediately lead to new top/bottom or sidewall scattering.
    This is necessary to prevent phonons leaving the structure boundaries."""
    x, y, z = move(ph, cf.timestep)
    return True if (abs(z) < cf.thickness / 2 and abs(x) < cf.width / 2 and y > 0) else False


def side_wall_scattering(ph, scattering_types):
    """Check if the phonon hits a side wall and output new vector"""
    x, y, z = move(ph, cf.timestep)

    # If phonon is beyond the side wall:
    if abs(x) > cf.width/2:

        # Calculate angle to the surface and specular scattering probability:
        a = acos(cos(ph.phi)*sin(abs(ph.theta))) # Angle to the surface
        p = specularity(a, cf.side_wall_roughness, ph.wavelength)

        # Specular scattering:
        if random() < p:
            scattering_types.walls = Scattering.SPECULAR
            ph.theta = - ph.theta

        # Diffuse scattering:
        else:
            scattering_types.walls = Scattering.DIFFUSE
            attempt = 0
            while attempt < 10:
                attempt += 1

                # Random distribution:
                # theta = -sign(x)*pi+sign(x)*pi*random()
                # phi = asin(2*random() - 1)

                # Lambert cosine distribution:
                ph.theta = -sign(x)*pi/2 + asin(2*random() - 1)
                ph.phi = asin((asin(2*random() - 1))/(pi/2))

                # Accept the angles if they do not cause new scattering:
                if no_new_scattering(ph):
                    break


def top_scattering(ph, scattering_types):
    """Check if the phonon hits the top surface and output new vector"""
    x, y, z = move(ph, cf.timestep)

    # If phonon is above the top surface, scattering happens:
    if z > cf.thickness/2:

        # Calculate angle to the surface and specular scattering probability:
        a = pi/2 - abs(ph.phi)
        p = specularity(a, cf.bottom_roughness, ph.wavelength)

        # Specular scattering:
        if random() < p:
            ph.phi = - ph.phi
            scattering_types.top_bottom = Scattering.SPECULAR

        # Diffuse scattering:
        else:
            # Random distribution:
            # theta = -pi + random()*2*pi
            # phi = -asin(random())

            # Lambert cosine distribution:
            ph.theta = - pi + random()*2*pi
            ph.phi = - random()*pi/2
            scattering_types.top_bottom = Scattering.DIFFUSE


def top_scattering_with_pillars(ph, scattering_types):
    """Check if the phonon hits the top surface and if this place has a pillar and output new vector"""
    x, y, z = move(ph, cf.timestep)

    # If phonon is below the bottom surface, scattering happens:
    if z > cf.thickness / 2:
        distances_from_centers = [0]*cf.pillar_coordinates.shape[0]
        for i in range(cf.pillar_coordinates.shape[0]):
            x0 = cf.pillar_coordinates[i,0]
            y0 = cf.pillar_coordinates[i,1]
            distances_from_centers[i] = (x-x0)**2 + (y-y0)**2

        # Angle to the surface:
        a = pi / 2 - abs(ph.phi)

        # If it is not under the pillar:
        if all((i > (cf.circular_hole_diameter / 2) ** 2) for i in distances_from_centers):
            p = specularity(a, cf.bottom_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.phi = - ph.phi
                scattering_types.top_bottom = Scattering.SPECULAR

            # Diffuse scattering:
            else:
                # Random distribution:
                # theta = -pi + random()*2*pi
                # phi = -asin(random())

                # Lambert cosine distribution:
                ph.theta = - pi + random()*2*pi
                ph.phi = - random()*pi/2
                scattering_types.top_bottom = Scattering.DIFFUSE

        # If it is the pillar top:
        elif z > cf.pillar_height + cf.thickness/2 and any((i > (cf.circular_hole_diameter / 2)**2) for i in distances_from_centers):
            p = specularity(a, cf.pillar_top_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.phi = - ph.phi
                scattering_types.top_bottom = Scattering.SPECULAR

            # Diffuse scattering:
            else:
                # Random distribution:
                # theta = -pi + random()*2*pi
                # phi = -asin(random())

                # Lambert cosine distribution:
                ph.theta = - pi + random()*2*pi
                ph.phi = - random()*pi/2
                scattering_types.top_bottom = Scattering.DIFFUSE


def bottom_scattering(ph, scattering_types):
    """Check if the phonon hits the bottom surface and calculate new angles"""
    x, y, z = move(ph, cf.timestep)

    # If phonon is below the top surface:
    if z < -cf.thickness/2:

        # Calculate angle to the surface and specular scattering probability:
        a = pi/2 - abs(ph.phi)
        p = specularity(a, cf.bottom_roughness, ph.wavelength)

        # Specular scattering:
        if random() < p:
            ph.phi = -ph.phi
            scattering_types.top_bottom = Scattering.SPECULAR

        # Diffuse scattering:
        else:
            # Random distribution:
            # theta = -pi + random()*2*pi
            # phi = asin(random())

            # Lambert cosine distribution:
            ph.theta = - pi + random()*2*pi
            ph.phi = random()*pi/2
            scattering_types.top_bottom = Scattering.DIFFUSE


def surface_scattering(ph, scattering_types):
    """Check if there will be a surface scattering on this timestep and return new direction"""

    # Scattering on top surface with and without pillars:
    if cf.include_pillars:
        top_scattering_with_pillars(ph, scattering_types)
    else:
        top_scattering(ph, scattering_types)

    # Scattering on bottom surface:
    if scattering_types.top_bottom is None:
        bottom_scattering(ph, scattering_types)

    # Scattering on sidewalls:
    side_wall_scattering(ph, scattering_types)

    # Scattering on parabolic walls:
    if cf.include_top_parabola:
        top_parabola_scattering(ph, scattering_types)
    if cf.include_bottom_parabola:
        bottom_parabola_scattering(ph, scattering_types)

    # Scattering on holes:
    if cf.include_holes:
        # Preliminary move to see if phonon would cross something:
        x, y, z = move(ph, cf.timestep)

        # Check for each hole:
        for i in range(cf.hole_coordinates.shape[0]):

            # Coordinates of the hole center:
            x0 = cf.hole_coordinates[i, 0]
            y0 = cf.hole_coordinates[i, 1]

            if cf.hole_shapes[i] == "circle":
                rad = cf.circular_hole_diameter * (1 + cf.hole_coordinates[i, 2]) / 2
                scattering_on_circular_holes(ph, x0, y0, rad, scattering_types, x, y, z)

            elif cf.hole_shapes[i] == "rectangle":
                # Correction of the hole size if there are holes of non-standard size:
                Lx = cf.rectangular_hole_side_x * (cf.hole_coordinates[i, 2] + 1)
                Ly = cf.rectangular_hole_side_y * (cf.hole_coordinates[i, 2] + 1)
                scattering_on_rectangular_holes(ph, x0, y0, Lx, Ly, scattering_types, x, y, z)

            elif cf.hole_shapes[i] == "triangle_up":
                Lx = cf.rectangular_hole_side_x * (cf.hole_coordinates[i, 2] + 1)
                Ly = cf.rectangular_hole_side_y * (cf.hole_coordinates[i, 2] + 1)
                scattering_on_triangle_up_holes(ph, x0, y0, Lx, Ly, scattering_types, x, y, z)

            elif cf.hole_shapes[i] == "triangle_down":
                Lx = cf.rectangular_hole_side_x * (cf.hole_coordinates[i, 2] + 1)
                Ly = cf.rectangular_hole_side_y * (cf.hole_coordinates[i, 2] + 1)
                scattering_on_triangle_down_holes(ph, x0, y0, Lx, Ly, scattering_types, x, y, z)

            # If there was any scattering, then no need to check other holes:
            if scattering_types.holes is not None:
                break

    # Scattering on pillars:
    if cf.include_pillars:

        # Preliminary move to see if phonon would cross something:
        x, y, z = move(ph, cf.timestep)

        for i in range(cf.pillar_coordinates.shape[0]):

            # Coordinates and radius of the given pillar:
            x0 = cf.pillar_coordinates[i, 0]
            y0 = cf.pillar_coordinates[i, 1]
            rad = cf.circular_hole_diameter * (1 + cf.pillar_coordinates[i,2]) / 2

            scattering_on_circular_pillars(ph, x0, y0, rad, scattering_types, x, y, z)

            # If there was any scattering, then no need to check other pillars:
            if scattering_types.pillars is not None:
                break

    # Correct angle if it became more than 180 degrees:
    ph.correct_angle()
