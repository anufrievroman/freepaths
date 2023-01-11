"""Modules provides scattering processes on various objects"""

from math import pi, cos, sin, tan, exp, sqrt, atan, asin, acos
from random import random
from numpy import sign

from move import move
from parameters import *
from lattices import hole_coordinates, pillar_coordinates, hole_shapes


def specularity(angle, roughness, wavelength):
    """Calculate probability of specular scattering with Soffer's equation"""
    return exp(-16*(pi**2)*(roughness**2)*((cos(angle))**2)/(wavelength**2))


def internal_scattering(ph, flight, scattering_types):
    """This function is checking if the time passed since previous diffuse scattering event
    reached the time until an internal scattering event, and if yes, scatters randomly"""
    if flight.time_since_previous_scattering >= ph.time_of_internal_scattering:
        ph.theta = -pi + random()*2*pi
        ph.phi = asin(2*random() - 1)
        scattering_types.internal = 'diffuse'


def reinitialization(ph, scattering_types):
    """Rethermalize phonon if it comes back to the hot side"""
    _, y, _ = move(ph)

    # If phonon returns to the staring line y = 0, generate it again:
    if y < 0:
        ph.assign_initial_coordinates()
        ph.assign_angles()
        scattering_types.hot_side = 'diffuse'


def scattering_on_circular_holes(ph, x0, y0, R, scattering_types):
    """Check if a phonon strikes a circular hole and calculate the new direction"""
    x, y, z = move(ph)

    # If phonon is inside the circle with radius R:
    if (x - x0)**2 + (y - y0)**2 <= R**2:

        # Calculate specular scattering probability (Soffer's equation):
        if y == y0: y += 1e-9 # Prevent division by zero
        tangent_theta = atan((x - x0)/(y - y0))
        a = acos(cos(ph.phi)*cos(ph.theta + sign(y - y0)*tangent_theta))  # Angle to the surface
        p = specularity(a, HOLE_ROUGHNESS, ph.wavelength)

        # Specular scattering:
        if random() < p:
            ph.theta = - ph.theta - pi + 2*tangent_theta
            scattering_types.holes = 'specular'

        # Diffuse scattering:
        else:
            scattering_types.holes = 'diffuse'
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


def scattering_on_rectangular_holes(ph, x0, y0, Lx, Ly, scattering_types):
    """Check if the phonon strikes a rectangular hole and calculate new direction"""
    x, y, z = move(ph)

    # If the phonon is inside the rectangle:
    if (abs(x - x0) <= Lx / 2) and (abs(y - y0) <= Ly / 2):

        # Coordinate y of the intersection with the hole side:
        y1 = (y0 - y) + cos(ph.theta)*(Lx/2 - abs(x0 - x))/abs(sin(ph.theta))

        # Scattering on left and right walls of the hole:
        if abs(y1) <= Ly/2:

            # Specular scattering probability (Soffer's equation):
            a = acos(cos(ph.phi)*sin(abs(ph.theta)))  # Angle to the normal to the surface
            p = specularity(a, HOLE_ROUGHNESS, ph.wavelength)

            # Specular scattering:
            if random() < p:
                scattering_types.holes = 'specular'
                ph.theta = - ph.theta

            # Diffuse scattering:
            else:
                scattering_types.holes = 'diffuse'
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
            # Specular scattering probability (Soffer's equation):
            a = acos(cos(ph.phi)*cos(ph.theta))     # Angle to the surface
            p = specularity(a, HOLE_ROUGHNESS, ph.wavelength)

            # Specular scattering:
            if random() < p:
                scattering_types.holes = 'specular'
                ph.theta = sign(ph.theta)*pi - ph.theta

            # Diffuse scattering:
            else:
                scattering_types.holes = 'diffuse'
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


def scattering_on_circular_pillars(ph, x0, y0, R_base, scattering_types):
    """Check if a phonon strikes a circular pillar and calculate new direction"""
    x, y, z = move(ph)

    # Cone radius at a given z coordinate:
    R = R_base - (z - THICKNESS / 2) / tan(PILLAR_WALL_ANGLE)

    # If phonon crosses the pillar boundary. Third condition is to exclude all other pillars:
    if ((x - x0)**2 + (y - y0)**2 >= R**2) and (z > THICKNESS / 2) \
            and ((x - x0) ** 2 + (y - y0) ** 2 < (R + 2 * ph.speed * TIMESTEP) ** 2):

        # Calculate specular scattering probability (Soffer's equation):
        tangent_theta = atan((x - x0)/(y - y0))
        a = atan(tan((pi/2 - ph.theta) + tangent_theta) * cos(ph.phi - (pi / 2 - PILLAR_WALL_ANGLE)))  # Angle to the surface
        p = specularity(a, PILLAR_ROUGHNESS, ph.wavelength)

        # Specular scattering:
        if random() < p:

            # If phonon moves from the center of the pillar to the wall:
            if sqrt((abs(x) - abs(x0))**2 + (abs(y) - abs(y0))**2) >= sqrt((abs(ph.x) - abs(x0))**2 + (abs(ph.y) - abs(y0))**2) :

                # If theta does not reflect back:
                if ph.phi < pi/2 - 2*PILLAR_WALL_ANGLE:
                    ph.phi = ph.phi - (pi / 2 - PILLAR_WALL_ANGLE)

                # Regular reflection:
                else:
                    ph.theta = - ph.theta - pi + 2*tangent_theta
                    ph.phi = ph.phi - (pi / 2 - PILLAR_WALL_ANGLE)

            # If phonon strikes the wall as it goes towards the center:
            else:
                ph.phi = -sign(ph.phi) * ph.phi - 2 * PILLAR_WALL_ANGLE
            scattering_types.pillars = 'specular'

        # Diffuse scattering:
        else:
            ph.theta = tangent_theta - asin(2*random()-1) + pi*(y >= y0)
            ph.phi = asin((asin(2*random() - 1))/(pi/2)) - (pi / 2 - PILLAR_WALL_ANGLE)
            scattering_types.pillars = 'diffuse'


def scattering_on_triangle_down_holes(ph, x0, y0, Lx, Ly, scattering_types):
    """Check if the phonon strikes a reverse triangular hole and calculate new direction after the scattering"""
    x, y, z = move(ph)

    # Angle of the triangle:
    beta = atan(0.5*Lx/Ly)

    # If phonon is inside the triangle:
    if (Ly/2 - (y - y0) <= (Lx/2 - abs(x - x0))/tan(beta)) and (abs(y - y0) < Ly/2):

        # Scattering on the bottom wall of the triangle:
        if (y + TIMESTEP * ph.speed > y0 + Ly / 2) and (abs(ph.theta) > pi / 2):

            # Calculate specular scattering probability (Soffer's equation):
            a = acos(cos(ph.phi)*cos(ph.theta)) # Angle to the surface
            p = specularity(a, HOLE_ROUGHNESS, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = sign(ph.theta)*pi - ph.theta
                scattering_types.holes = 'specular'

            # Diffuse scattering:
            else:
                # Lambert cosine distribution:
                ph.theta = asin(2*random() - 1)
                ph.phi = asin((asin(2*random() - 1))/(pi/2))
                scattering_types.holes = 'diffuse'

        # Scattering on the sidewalls of the triangle:
        else:
            # Calculate specular scattering probability (Soffer's equation):
            a = acos(cos(ph.phi)*cos(ph.theta - sign(x - x0)*(pi/2 - beta)))  # Angle to the surface
            p = specularity(a, HOLE_ROUGHNESS, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = - ph.theta + sign(x - x0)*2*beta
                scattering_types.holes = 'specular'

            # Diffuse scattering:
            else:
                rand_sign = sign((2*random() - 1))
                ph.theta = rand_sign*pi - rand_sign*asin(random()) - sign(x-x0)*(pi/2 - beta)
                ph.phi = asin((asin(2*random() - 1))/(pi/2))
                scattering_types.holes= 'diffuse'


def scattering_on_triangle_up_holes(ph, x0, y0, Lx, Ly, scattering_types):
    """Check if the phonon strikes a reverse triangular hole and calculate new direction"""
    x, y, z = move(ph)

    # Angle of the triangle:
    beta = atan(0.5*Lx/Ly)

    # If phonon is inside the triangle:
    if (Ly/2 + (y - y0) <= (Lx/2 - abs(x - x0))/tan(beta)) and (abs(y - y0) < Ly/2):
        #x1=Ly/2/tan(theta) - abs(y0-y)/tan(theta) + abs(x0-x)

        # Scattering on the bottom wall of the triangle:
        if ((y - TIMESTEP * ph.speed) < (y0 - Ly / 2)) and (abs(ph.theta) < pi / 2):

            # Calculate specular scattering probability (Soffer's equation):
            a = acos(cos(ph.phi)*cos(ph.theta)) # Angle to the surface
            p = specularity(a, HOLE_ROUGHNESS, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = sign(ph.theta)*pi - ph.theta
                scattering_types.holes = 'specular'

            # Diffuse scattering:
            else:
                # Lambert cosine distribution:
                rand_sign = sign((2*random() - 1))
                ph.theta = rand_sign*pi/2 + rand_sign*acos(random())
                ph.phi = asin((asin(2*random() - 1))/(pi/2))
                scattering_types.holes = 'diffuse'

        # Scattering on the sidewalls of the triangle:
        else:

            # Calculate specular scattering probability (Soffer's equation):
            a = acos(cos(ph.phi)*cos(ph.theta + sign(x - x0)*(pi/2 - beta)))  # Angle to the surface
            p = specularity(a, HOLE_ROUGHNESS, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = - ph.theta - sign(x - x0)*2*beta
                scattering_types.holes = 'specular'

            # Diffuse scattering:
            else:
                # Lambert cosine distribution:
                ph.theta = asin(2*random() - 1) + sign(x - x0)*(pi/2 - beta)
                ph.phi = asin((asin(2*random() - 1))/(pi/2))
                scattering_types.holes = 'diffuse'


def no_new_scattering(ph):
    """Check if new angles do not immediately lead to new top/bottom or sidewall scattering.
    This is necessary to prevent phonons leaving the structure boundaries."""
    x, y, z = move(ph)
    accept_angles = True if (abs(z) < THICKNESS / 2 and abs(x) < WIDTH / 2 and y > 0) else False
    return accept_angles


def side_wall_scattering(ph, scattering_types):
    """Check if the phonon hits a side wall and output new vector"""
    x, y, z = move(ph)

    # If phonon is beyond the side wall:
    if abs(x) > WIDTH/2:

        # Calculate specular scattering probability (Soffer's equation):
        a = acos(cos(ph.phi)*sin(abs(ph.theta))) # Angle to the surface
        p = specularity(a, SIDE_WALL_ROUGHNESS, ph.wavelength)

        # Specular scattering:
        if random() < p:
            scattering_types.walls = 'specular'
            ph.theta = - ph.theta

        # Diffuse scattering:
        else:
            scattering_types.walls = 'diffuse'
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
    x, y, z = move(ph)

    # If phonon is above the top surface, scattering happens:
    if z > THICKNESS/2:

        # Calculate specular scattering probability (Soffer's equation):
        a = pi/2 - abs(ph.phi)
        p = specularity(a, BOTTOM_ROUGHNESS, ph.wavelength)

        # Specular scattering:
        if random() < p:
            ph.phi = - ph.phi
            scattering_types.top_bottom = 'specular'

        # Diffuse scattering:
        else:
            # Random distribution:
            # theta = -pi + random()*2*pi
            # phi = -asin(random())

            # Lambert cosine distribution:
            ph.theta = - pi + random()*2*pi
            ph.phi = - random()*pi/2
            scattering_types.top_bottom = 'diffuse'


def top_scattering_with_pillars(ph, pillar_coordinates, scattering_types):
    """Check if the phonon hits the top surface and if this place has a pillar and output new vector"""
    x, y, z = move(ph)

    # If phonon is below the bottom surface, scattering happens:
    if z > THICKNESS/2:
        distances_from_centers = [0]*pillar_coordinates.shape[0]
        for i in range(pillar_coordinates.shape[0]):                            # For each pillar
            x0 = pillar_coordinates[i,0]                                          # Coordinates of the pillar center
            y0 = pillar_coordinates[i,1]
            distances_from_centers[i] = (x-x0)**2 + (y-y0)**2

        # Angle to the surface:
        a = pi / 2 - abs(ph.phi)

        # If it is not under the pillar:
        if all((i > (CIRCULAR_HOLE_DIAMETER / 2) ** 2) for i in distances_from_centers):
            p = specularity(a, BOTTOM_ROUGHNESS, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.phi = - ph.phi
                scattering_types.top_bottom = 'specular'

            # Diffuse scattering:
            else:
                # Random distribution:
                # theta = -pi + random()*2*pi
                # phi = -asin(random())

                # Lambert cosine distribution:
                ph.theta = - pi + random()*2*pi
                ph.phi = - random()*pi/2
                scattering_types.top_bottom = 'diffuse'

        # If it is the pillar top:
        elif z > PILLAR_HEIGHT + THICKNESS/2 and any((i > (CIRCULAR_HOLE_DIAMETER / 2)**2) for i in distances_from_centers):
            p = specularity(a, PILLAR_TOP_ROUGHNESS, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.phi = - ph.phi
                scattering_types.top_bottom = 'specular'

            # Diffuse scattering:
            else:
                # Random distribution:
                # theta = -pi + random()*2*pi
                # phi = -asin(random())

                # Lambert cosine distribution:
                ph.theta = - pi + random()*2*pi
                ph.phi = - random()*pi/2
                scattering_types.top_bottom = 'diffuse'


def bottom_scattering(ph, scattering_types):
    """Check if the phonon hits the bottom surface and calculate new angles"""
    x, y, z = move(ph)

    # If phonon is below the top surface:
    if z < -THICKNESS/2:

        # Calculate specular scattering probability (Soffer's equation):
        a = pi/2 - abs(ph.phi)
        p = specularity(a, BOTTOM_ROUGHNESS, ph.wavelength)

        # Specular scattering:
        if random() < p:
            ph.phi = -ph.phi
            scattering_types.top_bottom = 'specular'

        # Diffuse scattering:
        else:
            # Random distribution:
            # theta = -pi + random()*2*pi
            # phi = asin(random())

            # Lambert cosine distribution:
            ph.theta = - pi + random()*2*pi
            ph.phi = random()*pi/2
            scattering_types.top_bottom = 'diffuse'


def surface_scattering(ph, scattering_types):
    """This function checks if there will be a surface scattering on this timestep and returns a new direction"""

    # Scattering on top surface:
    if INCLUDE_PILLARS:
        top_scattering_with_pillars(ph, pillar_coordinates, scattering_types)
    else:
        top_scattering(ph, scattering_types)

    # Scattering on bottom surface:
    if scattering_types.top_bottom is None:
        bottom_scattering(ph, scattering_types)

    # Scattering on sidewalls:
    side_wall_scattering(ph, scattering_types)

    # Scattering on holes:
    if INCLUDE_HOLES:
        for i in range(hole_coordinates.shape[0]):

            # Coordinates of the hole center:
            x0 = hole_coordinates[i, 0]
            y0 = hole_coordinates[i, 1]

            if hole_shapes[i] == 'circle':
                rad = CIRCULAR_HOLE_DIAMETER * (1 + hole_coordinates[i, 2]) / 2
                scattering_on_circular_holes(ph, x0, y0, rad, scattering_types)

            elif hole_shapes[i] == 'rectangle':
                # Correction of the hole size if there are holes of non-standard size:
                Lx = RECTANGULAR_HOLE_SIDE_X * (hole_coordinates[i, 2] + 1)
                Ly = RECTANGULAR_HOLE_SIDE_Y * (hole_coordinates[i, 2] + 1)
                scattering_on_rectangular_holes(ph, x0, y0, Lx, Ly, scattering_types)

            elif hole_shapes[i] == 'triangle_down':
                Lx = RECTANGULAR_HOLE_SIDE_X * (hole_coordinates[i, 2] + 1)
                Ly = RECTANGULAR_HOLE_SIDE_Y * (hole_coordinates[i, 2] + 1)
                scattering_on_triangle_down_holes(ph, x0, y0, Lx, Ly, scattering_types)

            elif hole_shapes[i] == 'triangle_up':
                Lx = RECTANGULAR_HOLE_SIDE_X * (hole_coordinates[i, 2] + 1)
                Ly = RECTANGULAR_HOLE_SIDE_Y * (hole_coordinates[i, 2] + 1)
                scattering_on_triangle_up_holes(ph, x0, y0, Lx, Ly, scattering_types)

            # If there was any scattering, then no need to check other holes:
            if scattering_types.holes is not None:
                break

    # Scattering on pillars:
    if INCLUDE_PILLARS:
        for i in range(pillar_coordinates.shape[0]):

            # Coordinates and radius of the given pillar:
            x0 = pillar_coordinates[i, 0]
            y0 = pillar_coordinates[i, 1]
            rad = CIRCULAR_HOLE_DIAMETER * (1 + pillar_coordinates[i,2]) / 2

            scattering_on_circular_pillars(ph, x0, y0, rad, scattering_types)

            # If there was any scattering, then no need to check other pillars:
            if scattering_types.pillars is not None:
                break

    # Correct angle if it became more than 180 degrees:
    ph.correct_angle()
