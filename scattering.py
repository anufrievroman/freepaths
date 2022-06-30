import matplotlib.pyplot as plt
import numpy as np

from numpy import sign
from math import pi, cos, sin, tan, exp, log, sqrt, atan, asin, acos
from random import random, choice, randint
from scipy.constants import k, hbar
from functools import lru_cache

from main import initialization
from move import *
from parameters import *
from lattices import *


def scattering_on_rectangular_holes(x_orig, y_orig, z_orig, theta, phi, f, speed, x0, y0, Lx, Ly):
    '''This function checks if the phonon strikes a rectangular hole and what is the new direction after the scattering'''
    x, y, z = move(x_orig, y_orig, z_orig, theta, phi, speed)
    scattering_type = 'no_scattering'

    # If the phonon is inside the rectangle:
    if (abs(x - x0) <= Lx / 2) and (abs(y - y0) <= Ly / 2):

        # Coordinate y of the intersection with the hole side:
        y1 = (y0 - y) + cos(theta)*(Lx/2 - abs(x0 - x))/abs(sin(theta))

        # Scattering on left and right walls of the hole:
        if abs(y1) <= Ly/2:

            # Specular scattering probability (Soffer's equation):
            a = acos(cos(phi)*sin(abs(theta)))  # Angle to the normal to the surface
            p = exp(-16*(pi**2)*(hole_roughness**2)*((cos(a))**2)/((speed/f)**2))

            # Specular scattering:
            if random() < p:
                scattering_type = 'specular'
                theta = -theta

            # Diffuse scattering:
            else:
                scattering_type='diffuse'
                attempt = 0
                while attempt < 10:
                    attempt += 1

                    # Lambert cosine distribution:
                    theta = -sign(sin(theta))*pi/2 + asin(2*random()-1)
                    phi = asin((asin(2*random() - 1))/(pi/2))

                    # Accept the angles only if they do not lead to new scattering:
                    if no_new_scattering(x_orig, y_orig, z_orig, theta, phi, speed):
                        break

        # Scattering on top and bottom walls of the hole:
        else:
            # Specular scattering probability (Soffer's equation):
            a = acos(cos(phi)*cos(theta))                                       # Angle to the surface
            p = exp(-16*(pi**2)*(hole_roughness**2)*((cos(a))**2)/((speed/f)**2))

            # Specular scattering:
            if random() < p:
                scattering_type = 'specular'
                theta = sign(theta)*pi - theta

            # Diffuse scattering:
            else:
                scattering_type = 'diffuse'
                attempt = 0
                while attempt < 10:
                    attempt += 1

                    # Scattering on the top surface of the hole:
                    if abs(theta) > pi / 2:
                        # Lambert cosine distribution:
                        theta = asin(2*random() - 1)
                        phi = asin((asin(2*random() - 1))/(pi/2))

                    # Scattering on the bottom surface of the hole:
                    else:
                        # Lambert cosine distribution:
                        rand_sign = sign((2*random() - 1))
                        theta = rand_sign*pi/2 + rand_sign*acos(random())
                        phi = asin((asin(2*random() - 1))/(pi/2))

                    # Accept the angles only if they do not immediately cause new scattering:
                    if no_new_scattering(x_orig, y_orig, z_orig, theta, phi, speed):
                        break

    return theta, phi, scattering_type


def scattering_on_circular_holes(x_orig, y_orig, z_orig, theta, phi, f, speed, x0, y0, R):
    '''This function checks if a phonon strikes a circular hole and what is the new direction after the scattering'''
    x, y, z = move(x_orig, y_orig, z_orig, theta, phi, speed)
    scattering_type = 'no_scattering'

    # If phonon is inside the circle with radius R:
    if (x - x0)**2 + (y - y0)**2 <= R**2:

        # Calculate specular scattering probability (Soffer's equation):
        if y == y0: y += 1e-9 # Prevent division by zero
        tangent_theta = atan((x - x0)/(y - y0))
        a = acos(cos(phi)*cos(theta+sign(y-y0)*tangent_theta)) # Angle to the surface
        p = exp(-16*(pi**2)*(hole_roughness**2)*((cos(a))**2)/((speed/f)**2))

        # Specular scattering:
        if random() < p:
            scattering_type = 'specular'
            theta = -theta - pi + 2*tangent_theta
            phi = phi

        # Diffuse scattering:
        else:
            scattering_type='diffuse'
            attempt = 0
            while attempt < 10:
                attempt += 1

                # Random distribution:
                # theta = tangent_theta - (-pi/2 + pi*random()) - pi*(y < y0)
                # phi = asin(2*random() - 1)

                # Lambert cosine distribution:
                theta = tangent_theta - (asin(2*random() - 1)) - pi*(y < y0)
                phi = asin((asin(2*random() - 1))/(pi/2))

                # Accept the angles only if they do not immediately cause new scattering:
                if no_new_scattering(x_orig, y_orig, z_orig, theta, phi, speed):
                    break
    return theta, phi, scattering_type


def scattering_on_circular_pillars(x, y, z, theta, phi, frequency, speed, x0, y0, R_base):
    '''Check if a phonon strikes a circular pillar and calculate new direction'''
    x_previous = x
    y_previous = y
    x, y, z = move(x, y, z, theta, phi, speed)

    # Cone radius at a given z coordinate:
    R = R_base - (z - thickness/2)/tan(pillar_wall_angle)

    # Return phi into the -pi/2 : pi/2 range:
    phi = phi - sign(phi) * pi * (abs(phi) > pi / 2)

    # If phonon crosses the pillar boundary. Third condition is to exclude all other pillars:
    if ((x-x0)**2+(y-y0)**2 >= R**2) and (z > thickness/2) and ((x-x0)**2+(y-y0)**2 < (R+2*speed*timestep)**2):

        # Calculate specular scattering probability (Soffer's equation):
        tangent_theta = atan((x-x0)/(y-y0))
        lam = speed/frequency
        a = atan(tan((pi/2-theta)+tangent_theta)*cos(phi-(pi/2-pillar_wall_angle))) # Angle to the surface
        p = exp(-16*(pi**2)*(pillar_roughness**2)*((cos(pi/2-a))**2)/(lam**2))

        # Specular scattering:
        if random()<p:
            # If phonon moves from the center of the pillar to the wall:
            if sqrt((abs(x)-abs(x0))**2+(abs(y)-abs(y0))**2) >= sqrt((abs(x_previous)-abs(x0))**2+(abs(y_previous)-abs(y0))**2) :
                # If theta does not reflect back:
                if (phi < pi/2-2*pillar_wall_angle):
                    #new_phi=phi+2*pillar_wall_angle
                    new_phi=phi-(pi/2-pillar_wall_angle)
                    new_theta=theta
                # Regular reflection:
                else:
                    #new_phi=phi+2*pillar_wall_angle
                    new_theta=-theta-pi+2*tangent_theta
                    new_phi=phi-(pi/2-pillar_wall_angle)

            # If phonon strikes the wall as it goes towards the center:
            else:
                new_phi=-sign(phi)*phi-2*pillar_wall_angle
                new_theta=theta
            scattering_type='specular'

        # Diffuse scattering:
        else:
            new_theta = tangent_theta - asin(2*random()-1) + pi*(y >= y0)
            # new_phi=asin(2*random()-1)-(pi/2-pillar_wall_angle)
            new_phi = asin((asin(2*random() - 1))/(pi/2)) - (pi/2-pillar_wall_angle)
            scattering_type='diffuse'

    # No scattering
    else:
        new_theta = theta
        new_phi = phi
        scattering_type = 'no_scattering'
    return new_theta, new_phi, scattering_type


def scattering_on_triangle_down_holes(x_orig, y_orig, z_orig, theta, phi, f, speed, x0, y0, Lx, Ly):
    '''This function checks if the phonon strikes a reverse triangular hole and what is the new direction after the scattering'''
    x, y, z = move(x_orig, y_orig, z_orig, theta, phi, speed)
    scattering_type = 'no_scattering'

    # Angle of the triangle:
    beta = atan(0.5*Lx/Ly)

    # If phonon is inside the triangle:
    if (Ly/2 - (y - y0) <= (Lx/2 - abs(x - x0))/tan(beta)) and (abs(y - y0) < Ly/2):

        # Scattering on the bottom wall of the triangle:
        if (y + timestep*speed > y0 + Ly/2) and (abs(theta) > pi/2):

            # Calculate specular scattering probability (Soffer's equation):
            a = acos(cos(phi)*cos(theta)) # Angle to the surface
            p = exp(-16*(pi**2)*(hole_roughness**2)*((cos(a))**2)/((speed/f)**2))

            # Specular scattering:
            if random() < p:
                theta = sign(theta)*pi - theta
                phi = phi
                scattering_type = 'specular'

            # Diffuse scattering:
            else:
                # Lambert cosine distribution:
                theta = asin(2*random() - 1)
                phi = asin((asin(2*random() - 1))/(pi/2))
                scattering_type='diffuse'

        # Scattering on the sidewalls of the triangle:
        else:
            # Calculate specular scattering probability (Soffer's equation):
            a = acos(cos(phi)*cos(theta - sign(x - x0)*(pi/2 - beta)))  # Angle to the surface
            p = exp(-16*(pi**2)*(hole_roughness**2)*((cos(a))**2)/((speed/f)**2))

            # Specular scattering:
            if random() < p:
                theta = -theta + sign(x - x0)*2*beta
                phi = phi
                scattering_type = 'specular'

            # Diffuse scattering:
            else:
                rand_sign = sign((2*random() - 1))
                theta = rand_sign*pi - rand_sign*asin(random()) - sign(x-x0)*(pi/2 - beta)
                phi = asin((asin(2*random() - 1))/(pi/2))
                scattering_type='diffuse'

    return theta, phi, scattering_type


def scattering_on_triangle_up_holes(x_orig, y_orig, z_orig, theta, phi, f, speed, x0, y0, Lx, Ly):
    '''This function checks if the phonon strikes a reverse triangular hole and what is the new direction after the scattering'''
    x, y, z = move(x_orig, y_orig, z_orig, theta, phi, speed)
    scattering_type = 'no_scattering'

    # Angle of the triangle:
    beta = atan(0.5*Lx/Ly)

    # If phonon is inside the triangle:
    if (Ly/2 + (y - y0) <= (Lx/2 - abs(x - x0))/tan(beta)) and (abs(y - y0) < Ly/2):
        #x1=Ly/2/tan(theta) - abs(y0-y)/tan(theta) + abs(x0-x)

        # Scattering on the bottom wall of the triangle:
        if ((y-timestep*speed) < (y0-Ly/2)) and (abs(theta)<pi/2):

            # Calculate specular scattering probability (Soffer's equation):
            a = acos(cos(phi)*cos(theta)) # Angle to the surface
            p = exp(-16*(pi**2)*(hole_roughness**2)*((cos(a))**2)/((speed/f)**2))

            # Specular scattering:
            if random() < p:
                theta = sign(theta)*pi - theta
                phi = phi
                scattering_type = 'specular'

            # Diffuse scattering:
            else:
                # Lambert cosine distribution:
                rand_sign = sign((2*random() - 1))
                theta = rand_sign*pi/2 + rand_sign*acos(random())
                phi = asin((asin(2*random() - 1))/(pi/2))
                scattering_type='diffuse'

        # Scattering on the sidewalls of the triangle:
        else:

            # Calculate specular scattering probability (Soffer's equation):
            a = acos(cos(phi)*cos(theta+sign(x-x0)*(pi/2-beta)))  # Angle to the surface
            p = exp(-16*(pi**2)*(hole_roughness**2)*((cos(a))**2)/((speed/f)**2))

            # Specular scattering:
            if random()<p:
                scattering_type = 'specular'
                theta = -theta - sign(x - x0)*2*beta
                phi = phi

            # Diffuse scattering:
            else:
                # Lambert cosine distribution:
                theta = asin(2*random() - 1) + sign(x - x0)*(pi/2 - beta)
                phi = asin((asin(2*random() - 1))/(pi/2))
                scattering_type='diffuse'

    return theta, phi, scattering_type


def internal_scattering_time_calculation(frequency, polarization):
    '''Determine relaxation time after which this phonon will undergo internal scattering'''
    w = 2*pi*frequency                                              # Angular frequency

    # Depending on the material we assign different relaxation times:
    if material == "Si":
        deb_temp=152.0                                              # [K] Debay temperature
        tau_impurity=1/((2.95e-45)*(w**4))                          # Impurities scattering
        tau_umklapp=1/((0.95e-19)*(w**2)*T*exp(-deb_temp/T))        # Umklapp scattering
        #if polarization=='TA':
        #    tau_umklapp=4/((3.28e-19)*(w**2)*T*exp(-140/T))        # Ref. JAP 110, 034308 (2011)
        #if polarization=='LA':
        #    tau_umklapp=1/((3.28e-19)*(w**2)*T*exp(-140/T))
        #if polarization=='TA':
        #    tau_norm=1/(9.3e-13*w*(T**4))
        #if polarization=='LA':
        #    tau_norm=1/((2.0e-24)*(w**2)*(T**3))
        #tau_total=1/((1/tau_impurity)+(1/tau_norm)+(1/tau_umklapp))
        tau_internal = 1/((1/tau_impurity) + (1/tau_umklapp))

    if material == "SiC":
        # Parameters from Joshi et al, JAP 88, 265 (2000)
        deb_temp = 1200                                            # Debye temperature
        tau_impurity = 1/((8.46e-45)*(w**4.0))                     # Impurities scattering
        tau_umklapp = 1/((6.16e-20)*(w**2.0)*T*exp(-deb_temp/T))   # Umklapp scattering
        tau_4p = 1/((6.9e-23)*(T**2)*(w**2))                       # Four phonon processes
        tau_internal = 1/((1/tau_impurity) + (1/tau_umklapp) + (1/tau_4p))

    # Final relaxation time is determined with some randomization:
    time_of_internal_scattering = -log(random())*tau_internal      # Ref. PRB 94, 174303 (2016)
    return time_of_internal_scattering


def internal_scattering(theta, phi, time_since_previous_scattering, time_of_internal_scattering):
    '''This function is checking if the time passed since previous diffuse scattering event
    reached the time until an internal scattering event, and if yes, scatters randomly'''
    internal_scattering_type='none'
    if time_since_previous_scattering >= time_of_internal_scattering and internal_scattering_on:
        theta = -pi + random()*2*pi
        phi = asin(2*random()-1)
        internal_scattering_type = 'diffuse'
    return theta, phi, internal_scattering_type


def no_new_scattering(x, y, z, theta, phi, speed):
    '''Check if new angles do not immediately lead to new top/bottom or sidewall scattering.
    This is necessary to prevent phonons leaving the structure boundaries.'''
    x_new, y_new, z_new = move(x, y, z, theta, phi, speed)
    accept_angles = True if (abs(z_new) < thickness/2 and abs(x_new) < width/2 and y_new > 0) else False
    return accept_angles


def side_wall_scattering(x_orig, y_orig, z_orig, theta, phi, frequency, speed):
    '''Check if the phonon hits a side wall and output new vector'''
    x, y, z = move(x_orig, y_orig, z_orig, theta, phi, speed)
    scattering_type = 'no_scattering'

    # If phonon is beyond the side wall:
    if abs(x) > width/2:

        # Calculate specular scattering probability (Soffer's equation):
        lam = speed/frequency
        a = acos(cos(phi)*sin(abs(theta))) # Angle to the surface
        p = exp(-16*(pi**2)*(side_wall_roughness**2)*((cos(a))**2)/(lam**2))

        # Specular scattering:
        if random() < p:
            scattering_type = 'specular'
            theta = -theta

        # Diffuse scattering:
        else:
            scattering_type = 'diffuse'
            attempt = 0
            while attempt < 10:
                attempt += 1

                # Random distribution:
                # theta = -sign(x)*pi+sign(x)*pi*random()
                # phi = asin(2*random() - 1)

                # Lambert cosine distribution:
                theta = -sign(x)*pi/2 + asin(2*random() - 1)
                phi = asin((asin(2*random() - 1))/(pi/2))

                # Accept the angles only if they do not immediately cause new scattering:
                if no_new_scattering(x_orig, y_orig, z_orig, theta, phi, speed):
                    break
    return theta, phi, scattering_type


def top_scattering(x, y, z, theta, phi, frequency, speed):
    '''Check if the phonon hits the top surface and output new vector'''
    x, y, z = move(x, y, z, theta, phi, speed)
    scattering_type = 'no_scattering'

    # If phonon is above the top surface, scattering happens:
    if z > thickness/2:

        # Calculate specular scattering probability (Soffer's equation):
        lam = speed/frequency
        a = pi/2 - abs(phi)
        p = exp(-16*(pi**2)*(bottom_roughness**2)*((cos(a))**2)/(lam**2))

        # Specular scattering:
        if random() < p:
            phi = -phi
            scattering_type = 'specular'

        # Diffuse scattering:
        else:
            # Random distribution:
            # theta = -pi + random()*2*pi
            # phi = -asin(random())

            # Lambert cosine distribution:
            theta = -pi + random()*2*pi
            phi = -random()*pi/2
            scattering_type = 'diffuse'
    return theta, phi, scattering_type


def top_scattering_with_pillars(x, y, z, theta, phi, frequency, speed, pillar_coordinates):
    '''Check if the phonon hits the top surface and if this place has a pillar and output new vector'''
    x, y, z = move(x, y, z, theta, phi, speed)
    scattering_type = 'no_scattering'

    # If phonon is below the bottom surface, scattering happens:
    if z > thickness/2:
        distances_from_centers = [0]*pillar_coordinates.shape[0]
        for i in range(pillar_coordinates.shape[0]):                            # For each pillar
            x0 = pillar_coordinates[i,0]                                          # Coordinates of the pillar center
            y0 = pillar_coordinates[i,1]
            distances_from_centers[i] = (x-x0)**2 + (y-y0)**2

        # If it is not under the pillar:
        if all((i > (circular_hole_diameter/2)**2) for i in distances_from_centers):
            lam = speed/frequency
            a = pi/2 - abs(phi)
            p = exp(-16*(pi**2)*(bottom_roughness**2)*((cos(a))**2)/(lam**2))
            if random()<p:                                                      # Specular scattering
                phi = -phi
                scattering_type = 'specular'
            else:                                                               # Diffuse scattering
                # Random distribution:
                # theta = -pi + random()*2*pi
                # phi = -asin(random())

                # Lambert cosine distribution:
                theta = -pi + random()*2*pi
                phi = -random()*pi/2
                scattering_type = 'diffuse'

        # If it is the pillar top:
        elif z > pillar_height+thickness/2 and any((i > (circular_hole_diameter/2)**2) for i in distances_from_centers):
            lam=speed/frequency
            p=exp(-16*(pi**2)*(pillar_top_roughness**2)*((cos(pi/2-phi))**2)/(lam**2))     # Specular scattering probability
            if random()<p:                                                      # Specular scattering
                phi = -phi
                scattering_type = 'specular'
            else:                                                               # Diffuse scattering
                # Random distribution:
                # theta = -pi + random()*2*pi
                # phi = -asin(random())

                # Lambert cosine distribution:
                theta = -pi + random()*2*pi
                phi = -random()*pi/2
                scattering_type = 'diffuse'
    return theta, phi, scattering_type


def bottom_scattering(x, y, z, theta,phi, frequency, speed):
    '''Check if the phonon hits the bottom surface and calculate new angles'''
    x, y, z = move(x, y, z, theta, phi, speed)
    scattering_type = 'no_scattering'

    # If phonon is below the top surface:
    if z < -thickness/2:

        # Calculate specular scattering probability (Soffer's equation):
        lam = speed/frequency
        a = pi/2 - abs(phi)
        p = exp(-16*(pi**2)*(bottom_roughness**2)*((cos(a))**2)/(lam**2))

        # Specular scattering:
        if random() < p:
            phi = -phi
            scattering_type = 'specular'

        # Diffuse scattering:
        else:
            # Random distribution:
            # theta = -pi + random()*2*pi
            # phi = asin(random())

            # Lambert cosine distribution:
            theta = -pi + random()*2*pi
            phi = random()*pi/2
            scattering_type = 'diffuse'
    return theta, phi, scattering_type


def surface_scattering(x, y, z, theta, phi, frequency, hole_coordinates, hole_shapes, pillar_coordinates, speed):
    '''This function checks if there will be a surface scattering on this timestep and returns a new direction'''

    # Create variables assuming there is no scattering:
    hole_scattering_type = 'no_scattering'
    pillar_scattering_type = 'no_scattering'

    # Check for scattering on top surface:
    if pillars == 'yes':
        theta, phi, top_bottom_scattering_type = top_scattering_with_pillars(x, y, z, theta, phi, frequency, speed, pillar_coordinates)
    else:
        theta, phi, top_bottom_scattering_type = top_scattering(x, y, z, theta, phi, frequency, speed)

    # Check for scattering on bottom surface:
    if top_bottom_scattering_type == 'no_scattering':
        theta, phi, top_bottom_scattering_type = bottom_scattering(x, y, z, theta, phi, frequency, speed)

    # Check for scattering on sidewalls:
    theta, phi, wall_scattering_type = side_wall_scattering(x, y, z, theta, phi, frequency, speed)

    # Check for scattering on each hole:
    if holes == 'yes':
        for i in range(hole_coordinates.shape[0]):

            # Coordinates of the hole center:
            x0 = hole_coordinates[i, 0]
            y0 = hole_coordinates[i, 1]

            if hole_shapes[i] == 'circle':
                R = circular_hole_diameter*(1+hole_coordinates[i,2])/2
                # Scattering on this hole:
                theta,phi,hole_scattering_type=scattering_on_circular_holes(x,y,z,theta,phi,frequency,speed,x0,y0,R)

            elif hole_shapes[i] == 'rectangle':
                # Correction of the hole size if there are holes of non-standard size:
                Lx = rectangular_hole_side_x*(hole_coordinates[i,2]+1)
                Ly = rectangular_hole_side_y*(hole_coordinates[i,2]+1)
                # Scattering on this hole:
                theta,phi,hole_scattering_type=scattering_on_rectangular_holes(x,y,z,theta,phi,frequency,speed,x0,y0,Lx,Ly)

            elif hole_shapes[i] == 'triangle_down':
                Lx = rectangular_hole_side_x*(hole_coordinates[i,2]+1)
                Ly = rectangular_hole_side_y*(hole_coordinates[i,2]+1)
                theta,phi,hole_scattering_type=scattering_on_triangle_down_holes(x,y,z,theta,phi,frequency,speed,x0,y0,Lx,Ly)

            elif hole_shapes[i] == 'triangle_up':
                Lx = rectangular_hole_side_x*(hole_coordinates[i,2]+1)
                Ly = rectangular_hole_side_y*(hole_coordinates[i,2]+1)
                theta,phi,hole_scattering_type=scattering_on_triangle_up_holes(x,y,z,theta,phi,frequency,speed,x0,y0,Lx,Ly)

            # If there was any scattering, then no need to check other holes:
            if hole_scattering_type != 'no_scattering':
                break

    # Check for scattering on each pillar:
    if pillars == 'yes':
        for i in range(pillar_coordinates.shape[0]):

            # Coordinates and radius of the given pillar:
            x0 = pillar_coordinates[i,0]
            y0 = pillar_coordinates[i,1]
            R = circular_hole_diameter*(1 + pillar_coordinates[i,2])/2

            theta,phi,pillar_scattering_type = scattering_on_circular_pillars(x,y,z,theta,phi,frequency,speed,x0,y0,R)
            # If there was any scattering, then no need to check other pillars:
            if pillar_scattering_type != 'no_scattering':
                break

    # Check if angles are out of the [-pi:pi] range and return them back to this range:
    theta = theta - sign(theta)*2*pi*(abs(theta) > pi)
    # phi = phi - sign(phi)*2*pi*(abs(phi) > pi)

    # Finally, we record all scattering types that might have occurred and return it:
    surface_scattering_types = [wall_scattering_type, top_bottom_scattering_type, hole_scattering_type, pillar_scattering_type]
    return theta, phi, surface_scattering_types


def reinitialization(x, y, z, theta, phi, speed):
    '''Rethermalize phonon if it comes back to the hot side'''
    _, y_new, _ = move(x, y, z, theta, phi, speed)
    scattering_type = 'none'

    # If phonon returns to the staring line y = 0, generate it again:
    if y_new < 0:
        x, y, z, theta, phi = initialization()
        scattering_type = 'diffuse'
    return theta, phi, scattering_type, x, y, z

