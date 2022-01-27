from numpy import sign
from math import pi, cos, sin, tan, exp, log, sqrt, atan, asin, acos
import matplotlib.pyplot as plt
import numpy as np
from random import random, choice, randint
from scipy.constants import k, hbar
import os
import time
import sys
import shutil
from parameters import *
from lattices import hole_coordinates, hole_shapes, pillar_coordinates

plt.style.use('plotstyle.mplstyle')


def initialization():
    '''Assign initial position and angles of new phonons at the hot side'''

    # By default, phonons start at random position at the hot side:
    x = 0.49*width*(2*random()-1)
    y = 1e-12                       # It is practically zero
    z = 0.49*thickness*(2*random()-1)

    # In some lattices, phonons start from specific places:
    if hole_lattice_type=='serpentine':
        x = -width/2+(155e-9)/2 + 0.4*(155e-9)*(2*random()-1)
    elif hole_lattice_type=='diode_with_wires':
        x = 0.4*(period_x-rectangular_hole_side_x)*(2*random()-1)
    elif hole_lattice_type=='turn':
        x = 0.4*(period_x*5)*(2*random()-1)
    elif hole_lattice_type=='turn90':
        x = 0.4*width/2.0*(2*random()-1) - period_x*3.5
    elif hole_lattice_type=='directional_source':
        x = (2*random()-1)*50e-9
    elif hole_lattice_type=='tesla valve forward':
        x = (2*random()-1)*10e-9 - 100e-9
    elif hole_lattice_type=='tesla valve backward':
        x = (2*random()-1)*10e-9 - 100e-9

    # Choose distribution of initial angles: random, Lambert cosine, or directional
    if hot_side_angle_distribution == 'random':
        theta = -pi/2 + pi*random()
        phi = asin(2*random() - 1)
    if hot_side_angle_distribution == 'directional':
        theta = 0
        phi = -pi/2 + pi*random()
    if hot_side_angle_distribution == 'lambert':
        theta = asin(2*random() - 1)
        phi = asin((asin(2*random() - 1))/(pi/2))
    return x, y, z, theta, phi


def bulk_phonon_dispersion(N):
    '''Return phonon dispersion for given material calculated for N wavevectors over the G-X direction'''
    dispersion = np.zeros((N,4))
    if material == 'Si':                                            # Ref. APL 95 161901 (2009)
        dispersion[:,0] = [k*12e9/(N-1) for k in range(N)]                                                              # Wavevectors
        dispersion[:,1] = [abs(1369.42*k-2.405e-8*(k**2)-9.70e-19*(k**3)) for k in dispersion[:,0]]                     # LA branch
        dispersion[:,2] = [abs(1081.74*k-7.711e-8*(k**2)+5.674e-19*(k**3)+7.967e-29*(k**4)) for k in dispersion[:,0]]   # TA branch
        dispersion[:,3] = dispersion[:,2]                                                                               # TA branch
    if material == 'SiC':                                           # https://journals.aps.org/prb/pdf/10.1103/PhysRevB.50.17054
        dispersion[:,0] = [k*14414281503/(N-1) for k in range(N)]                                                       # Wavevectors
        dispersion[:,1] = [abs(-3.48834e-18*(k**3)+1.7604452e-08*(k**2)+1737.36296*k) for k in dispersion[:,0]]         # LA branch
        dispersion[:,2] = [abs(-2.21696e-19*(k**3)-3.4366886e-08*(k**2)+1077.98941*k) for k in dispersion[:,0]]         # TA branch
        dispersion[:,3] = dispersion[:,2]                                                                               # TA branch
    return dispersion


def phonon_properties_assignment():
    '''Assign phonon frequency (f) according to the Plank distribution at a given temperature T,
    chose polarization, and calculate group velocity from bulk dispersion'''
    if material == 'Si':
        default_speed = 6000                                          # [m/s] This is the speed for Debye approximation
    if material == 'SiC':
        default_speed = 6500                                          # [m/s] Need to change this probably!
    f_max=default_speed/(2*pi*hbar*default_speed/(2.82*k*T))          # Frequency of the peak of the Plank distribution
    DOS_max=3*((2*pi*f_max)**2)/(2*(pi**2)*(default_speed**3))        # DOS for f_max in Debye approximation
    bose_einstein_max=1/(exp((hbar*2*pi*f_max)/(k*T))-1)              # Bose-Einstein distribution for f_max
    plank_distribution_max=DOS_max*hbar*2*pi*f_max*bose_einstein_max  # Peak of the distribution (needed for normalization)

    # Load the phonon dispersion:
    dispersion = bulk_phonon_dispersion(1000)

    # Trying assign the frequency so that final frequency distribution would be Plankian:
    while True:                                                       # Until we obtain the frequency
        f=f_max*5*random()                                            # Let's draw a random frequency in the 0 - 5*f_max range
        DOS=3*((2*pi*f)**2)/(2*(pi**2)*(default_speed**3))            # Calculate the DOS in Debye approximation
        bose_einstein=1/(exp((hbar*2*pi*f)/(k*T))-1)                  # And the Bose-Einstein distribution
        plank_distribution=DOS*hbar*2*pi*f*bose_einstein              # Plank distribution

        # To statistically decide if this frequency belongs to our Plank distribution,
        # we take the normalized distribution at this frequency and draw a random number,
        # If the random number is lower than the distribution, we accept this frequency
        if random() < plank_distribution/plank_distribution_max and f < max(dispersion[:,1]):
            break

    # Depending on the branch, we calculate group velocity (speed) as dw/dk:
    polarization = choice(['TA','TA','LA'])                           # There are two TA branches and one LA branch
    if polarization == 'TA' and f < max(dispersion[:,2]):             # If TA polarization and frequency is good
        j=abs((np.abs(dispersion[:,2] - f)).argmin()-1)               # Calculate index closest to frequency f
        speed=2*pi*abs(dispersion[j+1,2]-dispersion[j,2])/abs(dispersion[j+1,0]-dispersion[j,0])
    else:                                                             # Otherwise it is LA polarization
        j=abs((np.abs(dispersion[:,1] - f)).argmin()-1)               # Calculate index closest to frequency f
        speed=2*pi*abs(dispersion[j+1,1]-dispersion[j,1])/abs(dispersion[j+1,0]-dispersion[j,0])
    phonon_properties = [f, polarization, speed]
    return phonon_properties


def phonon_properties_assignment_2(j, branch):
    '''This function assigns phonon frequency (f) according to the wavevector & branch and calculates group velocity from bulk dispersion'''
    dispersion=bulk_phonon_dispersion(number_of_phonons+1)
    K=(dispersion[j+1,0]+dispersion[j,0])/2                          # Wavevector (we take average in the interval)
    dK=(dispersion[j+1,0]-dispersion[j,0])                           # Delta wavevector
    w=2*pi*abs((dispersion[j+1,branch+1]+dispersion[j,branch+1])/2)  # Angular frequency  (we take average in the interval)
    dw=2*pi*abs(dispersion[j+1,branch+1]-dispersion[j,branch+1])     # Delta angular frequency
    speed = dw/dK                                                    # Group velocity
    frequency = w/(2*pi)
    polarization='LA' if branch == 0 else 'TA'
    phonon_properties = [frequency, polarization, speed]
    return phonon_properties, w, K, dK


# def heat_capacity_calculation():
#     '''This function calculated heat capacity from the phonon dispersion'''
#     specific_heat_capacity=0
#     for branch in range(3):                                                            # For each phonon branch
#         for j in range(number_of_phonons):
#             phonon_properties, w, K, dK  = phonon_properties_assignment_2(j,branch)
#             frequency, polarization, speed = phonon_properties
#             heat_capacity=k*((hbar*w/(k*T))**2)*exp(hbar*w/(k*T))/((exp(hbar*w/(k*T))-1)**2)    # Ref. PRB 88 155318 (2013)
#             specific_heat_capacity += heat_capacity
#     return specific_heat_capacity


def move(x, y, z, theta, phi, speed):
    '''This function moves a phonon in one timestep and returns new coordinates'''
    x += sin(theta)*abs(cos(phi))*speed*timestep
    y += cos(theta)*abs(cos(phi))*speed*timestep
    z += sin(phi)*speed*timestep
    return x, y, z


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


def distribution_calculation(filename, data_range, number_of_nodes):
    '''This function calculates distribution of numbers (histogram) in a given file'''
    data = np.loadtxt(filename)
    if data_range == None:
        data_range = np.max(data)
    distribution = np.zeros((number_of_nodes, 2))
    distribution[:,0] = np.linspace(0, data_range, number_of_nodes)
    distribution[:,1], _ = np.histogram(data[data != 0], number_of_nodes, range=(0, data_range))
    return distribution


def angle_distribution_calculation():
    '''This function analyses measured phonon angles and creates their distribution'''
    exit_angles = np.loadtxt("All exit angles.csv", dtype='float')
    initial_angles = np.loadtxt("All initial angles.csv", dtype='float')
    distribution=np.zeros((180,3))
    distribution[:,0]=range(-90,90)
    distribution[:,1], _ = np.histogram(np.degrees(exit_angles[exit_angles != 0]), 180, range=(-90, 90))
    distribution[:,2], _ = np.histogram(np.degrees(initial_angles), 180, range=(-90, 90))
    return distribution


def wavelength_distribution_calculation(number_of_nodes):
    '''This function calculates phonon wavelength distribution from their frequencies and velocities'''
    frequencies = np.loadtxt("All initial frequencies.csv")
    speeds = np.loadtxt("All group velocities.csv")
    wavelengths = np.zeros((len(speeds)))
    wavelengths[:] = speeds[:]/frequencies[:]
    data_range = np.amax(wavelengths)
    distribution = np.zeros((number_of_nodes, 2))
    distribution[:,0] = np.linspace(0, data_range, number_of_nodes)
    distribution[:,1], _ = np.histogram(wavelengths, number_of_nodes, range=(0, data_range))
    return distribution


def scattering_events_statistics_calculation(statistics_of_scattering_events, surface_scattering_types,
                                             reinitialization_scattering_type, internal_scattering_type, y):
    '''This function analyzes type of scattering events at the current timestep and adds them to the statistics'''
    # Calculate in which length segment (starting from zero) we are:
    segment = int(y//(length/number_of_length_segments))

    # Scattering on side walls:
    statistics_of_scattering_events[segment, 0] += 1 if surface_scattering_types[0] == 'diffuse' else 0
    statistics_of_scattering_events[segment, 1] += 1 if surface_scattering_types[0] == 'specular' else 0
    # Scattering on top and bottom:
    statistics_of_scattering_events[segment, 2] += 1 if surface_scattering_types[1] == 'diffuse' else 0
    statistics_of_scattering_events[segment, 3] += 1 if surface_scattering_types[1] == 'specular' else 0
    # Scattering on holes:
    statistics_of_scattering_events[segment, 4] += 1 if surface_scattering_types[2] == 'diffuse' else 0
    statistics_of_scattering_events[segment, 5] += 1 if surface_scattering_types[2] == 'specular' else 0
    # Internal scattering and rethermalization on hot side:
    statistics_of_scattering_events[segment, 6] += 1 if reinitialization_scattering_type == 'diffuse' else 0
    statistics_of_scattering_events[segment, 7] += 1 if internal_scattering_type == 'diffuse' else 0
    # Scattering on pillars
    statistics_of_scattering_events[segment, 8] += 1 if surface_scattering_types[3] == 'diffuse' else 0
    statistics_of_scattering_events[segment, 9] += 1 if surface_scattering_types[3] == 'specular' else 0
    return statistics_of_scattering_events


def scattering_map_calculation(x,y,scattering_maps,internal_scattering_type,surface_scattering_types):
    '''This function records the place where a scattering event occurred according to the event type'''
    if any(scattering == 'diffuse' for scattering in surface_scattering_types):
        scattering_maps[0].append(x)
        scattering_maps[1].append(y)
    elif any(scattering == 'specular' for scattering in surface_scattering_types):
        scattering_maps[2].append(x)
        scattering_maps[3].append(y)
    elif internal_scattering_type == 'diffuse':
        scattering_maps[4].append(x)
        scattering_maps[5].append(y)
    return scattering_maps


def record_time_in_segments(time_in_segments, y, group):
    '''Record how long phonons stay in different segments'''
    for segment_number in range(number_of_length_segments):
        segment_begining = segment_number*(length/number_of_length_segments)
        segment_end = (segment_number+1)*(length/number_of_length_segments)
        for ph in range(number_of_phonons_in_a_group):
            times = [timestep for point in y[:, ph] if ((segment_begining <= point < segment_end) and point != 0.0)]
            time_in_segments[segment_number, group] += sum(times)
    return time_in_segments


def maps_and_profiles_calculation(x,y,maps_and_profiles,phonon_properties,timestep_number,theta,phi):
    '''This function registers the phonon in the pixel corresponding to its current position and at certain timesteps
    and calculates thermal maps and thermal profiles along different axes'''
    # Unpacking stuff from the tuple container:
    thermal_map, heat_flux_profile_x, heat_flux_profile_y, temperature_profile_x, temperature_profile_y = maps_and_profiles
    frequency, polarization, speed = phonon_properties[0:3]

    # Calculate the index of the pixel in which this phonon is now:
    index_x = int(((x+width/2)*number_of_pixels_x) // width)
    # index_y = int((y*number_of_pixels_y) // length)
    index_y = int(y//(length/number_of_pixels_y))

    # Calculate the volume of this pixel:
    Vcell_x = length*thickness*width/number_of_pixels_x
    Vcell_y = width*thickness*length/number_of_pixels_y

    # Assign density depending on the material:
    if material == 'Si':
        material_density = 2330 # [kg/m^3]
    elif material == 'SiC':
        material_density = 3215 # [kg/m^3]

    # Here we arbitraraly correct the volume of the unit cells:
    if pillars=='yes':
        Vcell_x+=2.5*0.3333*pillar_height*(circular_hole_diameter/2)**2
        Vcell_y+=2.5*0.3333*pillar_height*(circular_hole_diameter/2)**2

    # Prevent error if the phonon is outside the structure:
    if (0 <= index_x < number_of_pixels_x) and (0 <= index_y < number_of_pixels_y):
        # Record energy h*w of this phonon into the pixel of thermal map:
        thermal_map[index_y,index_x] += hbar*2*pi*frequency

        # Record energy of this phonon into flux and temperature profiles:
        timeframe_number = int(((timestep_number+randint(0, number_of_timesteps))*timestep*number_of_timeframes) // (number_of_timesteps*timestep))
        if timeframe_number < number_of_timeframes:
            heat_flux_profile_x[index_x,timeframe_number]   += hbar*2*pi*frequency*cos(theta)*abs(cos(phi))*speed/Vcell_x
            heat_flux_profile_y[index_y,timeframe_number]   += hbar*2*pi*frequency*cos(theta)*abs(cos(phi))*speed/Vcell_y
            temperature_profile_x[index_x,timeframe_number] += hbar*2*pi*frequency/(specific_heat_capacity*material_density)/Vcell_x
            temperature_profile_y[index_y,timeframe_number] += hbar*2*pi*frequency/(specific_heat_capacity*material_density)/Vcell_y

    # Pack everything back into the tuple container and return:
    maps_and_profiles = [thermal_map, heat_flux_profile_x, heat_flux_profile_y, temperature_profile_x, temperature_profile_y]
    return maps_and_profiles


def create_empty_maps():
    '''This function just creates empty maps'''
    # Create empty profile map arrays:
    thermal_map=np.zeros((number_of_pixels_y,number_of_pixels_x))
    heat_flux_profile_x=np.zeros((number_of_pixels_x,number_of_timeframes))
    heat_flux_profile_y=np.zeros((number_of_pixels_y,number_of_timeframes))
    temperature_profile_x=np.zeros((number_of_pixels_x,number_of_timeframes))
    temperature_profile_y=np.zeros((number_of_pixels_y,number_of_timeframes))
    maps_and_profiles=[thermal_map,heat_flux_profile_x,heat_flux_profile_y,temperature_profile_x,temperature_profile_y]

    # Create empty scattering maps:
    diffuse_scattering_map_x  = []
    diffuse_scattering_map_y  = []
    specular_scattering_map_x = []
    specular_scattering_map_y = []
    internal_scattering_map_x = []
    internal_scattering_map_y = []

    # Package everything into one tuple:
    scattering_maps = [diffuse_scattering_map_x, diffuse_scattering_map_y, specular_scattering_map_x,
            specular_scattering_map_y, internal_scattering_map_x, internal_scattering_map_y]
    return maps_and_profiles, scattering_maps


def phonon_is_in_system(x, y):
    '''This function checks if the phonon at this timestep is still in the system and did not reach the cold side'''
    small_offset = 10e-9
    # Depending on where we set the cold side, we check if phonon crossed that line:
    if cold_side == 'top':
        phonon_is_in_system = (y < length)
    elif cold_side == 'right':
        phonon_is_in_system = ((y < length-1.1e-6) or (y > length-1.1e-6 and x < width/2.0-small_offset))
    elif cold_side == 'top and right':
        phonon_is_in_system = ((y < 1.0e-6) or (y > 1.0e-6 and x < width/2.0-small_offset and y < length))
    return phonon_is_in_system


def progress_bar(i, j, old_progress, scheme):
    '''This is a progress bar that outputs the progress each one percent but not more often'''
    if scheme == 1:
        progress = 100*(i*number_of_phonons_in_a_group+j)//number_of_phonons
    elif scheme == 2:
        progress = 100*(i*number_of_phonons+j)//(number_of_phonons*3)

    # Only update the progress bar if it changed by 1%
    if progress > old_progress:
        sys.stdout.write('\r'+'Progress: '+str(progress)+'%')
        sys.stdout.flush()
    return progress


def write_files(free_paths, free_paths_along_y, frequencies, exit_angles, initial_angles, group_velocities,
                statistics_of_scattering_events, all_travel_times, all_detected_frequencies):
    '''This function analyzes writes files with statistics'''
    sys.stdout.write('\r'+'Progress: 100%')
    sys.stdout.write("\n")

    os.chdir(output_folder_name)

    # Write the files with raw data:
    header = "Sidewalls diffuse, Sidewalls specular, Top & bottom diffuse, Top & bottom specular, Holes diffuse, Holes specular, Hot side, Internal, Pillars diffuse, Pillars specular"
    np.savetxt("Scattering events statistics.csv", statistics_of_scattering_events, fmt='%1.3e', delimiter=",", header=header)
    np.savetxt("All free paths.csv", free_paths, fmt='%1.3e', delimiter=",", header="L [m]")
    np.savetxt("All free paths in plane.csv", free_paths_along_y, fmt='%1.3e', delimiter=",", header="Ly [m]")
    np.savetxt("All initial frequencies.csv", frequencies, fmt='%1.3e', delimiter=",", header="f [Hz]")
    np.savetxt("All detected frequencies.csv", all_detected_frequencies, fmt='%1.3e', delimiter=",", header="f [Hz]")
    np.savetxt("All exit angles.csv", exit_angles, fmt='%1.3e', delimiter=",", header="Angle [rad]")
    np.savetxt("All initial angles.csv", initial_angles, fmt='%1.3e', delimiter=",", header="Angle [rad]")
    np.savetxt("All group velocities.csv", group_velocities, fmt='%1.3e', delimiter=",", header="Vg [rad]")
    np.savetxt("All travel times.csv", all_travel_times, fmt='%1.3e', delimiter=",", header="Travel time [s]")

    # Here we temporary save the united statistics on scattering in the entire structure
    united_statistics = []
    for i in range(10):
        united_statistics.append(sum([statistics_of_scattering_events[segment, i] for segment in range(number_of_length_segments)]))
    np.savetxt("Statistics.csv", united_statistics, fmt='%2.3e', delimiter=",")


def output_trajectories(x, y, z, N):
    '''This function outputs the phonon trajectories of N phonons'''
    # Create XY plot:
    fig, ax = plt.subplots()
    for i in range(N):
        ax.plot (np.trim_zeros(x[:,i])*1e6,np.trim_zeros(y[:,i])*1e6, linewidth=0.2)
    ax.set_xlabel('X (μm)', fontsize=12)
    ax.set_ylabel('Y (μm)', fontsize=12)
    ax.set_aspect('equal', 'datalim')
    fig.savefig("Phonon paths XY.pdf", dpi=300, format = 'pdf', bbox_inches="tight")
    if output_in_terminal: plt.show()

    # Create YZ plots:
    fig, ax = plt.subplots()
    for i in range(N):
        ax.plot (np.trim_zeros(y[:,i])*1e6,np.trim_zeros(z[:,i])*1e6, linewidth=0.1)
    ax.set_xlabel('Y (μm)', fontsize=12)
    ax.set_ylabel('Z (μm)', fontsize=12)
    ax.set_aspect('equal', 'datalim')
    fig.savefig("Phonon paths YZ.pdf", dpi=300, format = 'pdf', bbox_inches="tight")
    if output_in_terminal: plt.show()


def output_thermal_map(thermal_map):
    '''This function outputs the thermal map'''
    from matplotlib.colors import LogNorm
    minimum_of_colorbar=1e-20                                     # Cannot be zero!
    # Writing data into the file, if it was set in parameters:
    if output_raw_thermal_map:
        np.savetxt("Thermal map.csv", thermal_map, fmt='%1.2e', delimiter=",")
    thermal_map=np.flipud(thermal_map)
    #plt.imshow(thermal_map, cmap='hot', interpolation='none', extent=[(-width/2)*1e6,(width/2)*1e6,0,length*1e6] )                # can also use interpolation='bicubic' and norm=LogNorm(vmin=0.01, vmax=np.amax(thermal_map))
    fig = plt.figure()
    plt.imshow(thermal_map, cmap='hot', interpolation='none', extent=[(-width/2)*1e6,(width/2)*1e6,0,length*1e6], norm=LogNorm(vmin=minimum_of_colorbar, vmax=np.amax(thermal_map)) )
    plt.xlabel('X (μm)', fontsize=12)
    plt.ylabel('Y (μm)', fontsize=12)
    cbar=plt.colorbar()
    cbar.set_label('Energy density', rotation=90)
    fig.savefig("Thermal map.pdf", bbox_inches="tight")
    if output_in_terminal: plt.show()


def output_scattering_maps(scattering_maps):
    '''This function outputs scattering map of diffusive, specular, and internal scattering events'''
    fig, ax = plt.subplots()
    ax.plot (scattering_maps[2][:], scattering_maps[3][:], 'o', color='g', markersize=0.2, alpha=0.2)
    ax.plot (scattering_maps[4][:], scattering_maps[5][:], 'o', color='r', markersize=0.2, alpha=0.2)
    ax.plot (scattering_maps[0][:], scattering_maps[1][:], 'o', color='b', markersize=0.2, alpha=0.2)
    ax.set_xlabel('X (μm)', fontsize=12)
    ax.set_ylabel('Y (μm)', fontsize=12)
    ax.set_aspect('equal', 'datalim')
    fig.savefig("Scattering map.pdf", bbox_inches="tight")
    if output_in_terminal: plt.show()

    N = max(len(scattering_maps[i]) for i in range(6))   # Maximal number of scattering event of any type
    data = np.zeros((N,6))                               # Create a numpy array that we will output into a file
    for i in range(6):
        for value in scattering_maps[i]:
            data[j,i] = value
    np.savetxt("Scattering map.csv", data, fmt='%1.2e', delimiter=",")


def output_thermal_conductivity(maps_and_profiles):
    '''Calculate the thermal conductivity for each time interval from heat flux
    and temperature profiles accumulated in that interval'''
    thermal_map, J_profiles_x, J_profiles_y, T_profiles_x, T_profiles_y = maps_and_profiles
    thermal_conductivity = np.zeros((number_of_timeframes, 2))
    thermal_conductivity[:, 0] = range(number_of_timeframes)
    thermal_conductivity[:, 0] *= number_of_timesteps * timestep / number_of_timeframes

    # For each time interval calculate the thermal conductivity:
    for i in range(number_of_timeframes):
        # Here we ignore the first pixel, because there are anomalies usually...
        dT = T_profiles_y[1, i] - T_profiles_y[(number_of_pixels_y - 1), i]
        J = sum(J_profiles_y[1:number_of_pixels_y, i]) / (number_of_pixels_y - 1)
        # Here dL is shorter then aclual length because we ignore 1st pixel and lose one more due to averaging:
        dL = (number_of_pixels_y - 2) * length / number_of_pixels_y
        # By definition, J = -K*grad(T), so the thermal conductivity:
        thermal_conductivity[i, 1] = J * dL / dT

    # Create the plot
    fig, ax = plt.subplots()
    ax.plot(thermal_conductivity[:, 0] * 1e9, thermal_conductivity[:, 1], linewidth=1)
    ax.set_ylabel('Thermal conductivity (W/mK)', fontsize=12)
    ax.set_xlabel('Time (ns)', fontsize=12)
    fig.savefig("Thermal conductivity.pdf", dpi=300, format='pdf', bbox_inches="tight")
    if output_in_terminal: plt.show()
    np.savetxt('Thermal conductivity.csv', thermal_conductivity, fmt='%1.3e', delimiter=",")


def output_profiles(maps_and_profiles):
    '''Output temperature and flux profiles along the structure'''
    thermal_map, J_profiles_x, J_profiles_y, T_profiles_x, T_profiles_y = maps_and_profiles

    coordinates_x=np.arange(T_profiles_x.shape[0])*1e6*width/T_profiles_x.shape[0]      # Let's create coordinate arrays (in um)
    coordinates_y=np.arange(T_profiles_y.shape[0])*1e6*length/T_profiles_y.shape[0]

    # Saving all the profiles in the files:
    np.savetxt("Temperature profiles x.csv", np.vstack((coordinates_x,T_profiles_x.T)).T, fmt='%1.3e', delimiter=",")
    np.savetxt("Temperature profiles y.csv", np.vstack((coordinates_y,T_profiles_y.T)).T, fmt='%1.3e', delimiter=",")
    np.savetxt("Heat flux profiles x.csv", np.vstack((coordinates_x,J_profiles_x.T)).T, fmt='%1.3e', delimiter=",")
    np.savetxt("Heat flux profiles y.csv", np.vstack((coordinates_y,J_profiles_y.T)).T, fmt='%1.3e', delimiter=",")

    # Creating the plot of temperature profile:
    fig, ax = plt.subplots()
    for i in range(number_of_timeframes):
        ax.plot (coordinates_y[1:number_of_pixels_y],T_profiles_y[1:number_of_pixels_y,i], linewidth=1)
    ax.set_xlabel('Y (μm)', fontsize=12)
    ax.set_ylabel('Temperature (K)', fontsize=12)
    fig.savefig("Temperature profile.pdf", dpi=300, format = 'pdf', bbox_inches="tight")
    if output_in_terminal: plt.show()

    # Creating the plot of heat flux profile:
    fig, ax = plt.subplots()
    for i in range(number_of_timeframes):
        ax.plot (coordinates_y[0:number_of_pixels_y], J_profiles_y[0:number_of_pixels_y, i], linewidth=1)
    ax.set_xlabel('Y (μm)', fontsize=12)
    ax.set_ylabel('Heat flux (W/m^2)', fontsize=12)
    fig.savefig("Heat flux profile.pdf", dpi=300, format = 'pdf', bbox_inches="tight")
    if output_in_terminal: plt.show()


def output_distributions():
    '''This function outputs the distributions into the terminal and the folder'''
    angle_distributions = angle_distribution_calculation()
    frequency_distribution = distribution_calculation("All initial frequencies.csv", None, number_of_nodes)
    detected_frequency_distribution = distribution_calculation("All detected frequencies.csv", None, number_of_nodes)
    wavelength_distribution = wavelength_distribution_calculation(number_of_nodes)

    fig, ax = plt.subplots()
    ax.plot (angle_distributions[:,0],angle_distributions[:,1],'b')
    ax.plot (angle_distributions[:,0],angle_distributions[:,2],'r')
    ax.set_xlabel('Angle (degree)', fontsize=12)
    ax.set_ylabel('Number of phonons', fontsize=12)
    fig.savefig("Distribution of angles.pdf", dpi=300, format='pdf', bbox_inches="tight")
    if output_in_terminal: plt.show()
    np.savetxt('Distribution of angles.csv', angle_distributions, fmt='%1.3e', delimiter=",")

    free_path_distribution = distribution_calculation("All free paths.csv", length, number_of_nodes)
    fig, ax = plt.subplots()
    ax.plot (free_path_distribution[:,0]*1e6,free_path_distribution[:,1])
    ax.set_xlabel('Free flights (μm)', fontsize = 12)
    ax.set_ylabel('Number of flights', fontsize=12)
    fig.savefig("Distribution of free paths.pdf", dpi=300, format = 'pdf', bbox_inches="tight")
    if output_in_terminal: plt.show()
    np.savetxt('Distribution of free paths.csv', free_path_distribution, fmt='%1.3e', delimiter=",")

    fig, ax = plt.subplots()
    ax.plot (frequency_distribution[:,0],frequency_distribution[:,1])
    ax.set_xlabel('Frequency (Hz)', fontsize=12)
    ax.set_ylabel('Number of phonons', fontsize=12)
    fig.savefig("Distribution of initial frequencies.pdf", dpi=300, format = 'pdf', bbox_inches="tight")
    if output_in_terminal: plt.show()
    np.savetxt('Distribution of initial frequencies.csv', frequency_distribution, fmt='%1.3e', delimiter=",")

    fig, ax = plt.subplots()
    ax.plot (wavelength_distribution[:,0]*1e9,wavelength_distribution[:,1])
    ax.set_xlabel('Wavelength (nm)', fontsize=12)
    ax.set_ylabel('Number of phonons', fontsize=12)
    fig.savefig("Distribution of wavelengths.pdf", dpi=300, format = 'pdf', bbox_inches="tight")
    if output_in_terminal: plt.show()
    np.savetxt('Distribution of wavelengths.csv', wavelength_distribution, fmt='%1.3e', delimiter=",")

    travel_time_distribution = distribution_calculation("All travel times.csv", None, number_of_nodes)
    fig, ax = plt.subplots()
    ax.plot (travel_time_distribution[:,0]*1e9,travel_time_distribution[:,1])
    ax.set_xlabel('Travel time (ns)', fontsize=12)
    ax.set_ylabel('Number of phonons', fontsize=12)
    fig.savefig("Distribution of travel times.pdf", dpi=300, format = 'pdf', bbox_inches="tight")
    if output_in_terminal: plt.show()
    np.savetxt('Distribution of travel times.csv', travel_time_distribution, fmt='%1.3e', delimiter=",")

    fig, ax = plt.subplots()
    ax.plot (detected_frequency_distribution[:,0],detected_frequency_distribution[:,1])
    ax.set_xlabel('Frequency (Hz)', fontsize=12)
    ax.set_ylabel('Number of phonons', fontsize=12)
    fig.savefig("Distribution of detected frequencies.pdf", dpi=300, format = 'pdf', bbox_inches="tight")
    if output_in_terminal: plt.show()
    np.savetxt('Distribution of detected frequencies.csv', detected_frequency_distribution, fmt='%1.3e', delimiter=",")

    speeds = np.loadtxt("All group velocities.csv")
    frequencies = np.loadtxt("All initial frequencies.csv")
    fig, ax = plt.subplots()
    ax.plot (frequencies,speeds,'.')
    ax.set_xlabel('Frequency (Hz)', fontsize=12)
    ax.set_ylabel('Group velocity (m/s)', fontsize=12)
    fig.savefig("Group velocities.pdf", dpi=300, format = 'pdf', bbox_inches="tight")
    if output_in_terminal: plt.show()


def output_scattering_statistics(statistics_of_scattering_events):
    '''Calculate and plot rates of different scattering events in each length segment'''
    # Load the data from files:
    scattering_data = np.genfromtxt("Scattering events statistics.csv", unpack = True,  delimiter=',', skip_header = 1)
    segments, time_spent = np.genfromtxt("Time spent in segments.csv", unpack = True,  delimiter=',', skip_header = 1)

    # Create the plot:
    fig, ax = plt.subplots()
    all_scattering_rates = []
    for scattering_type in range(scattering_data.shape[0]):
        # Calculate scattering events per second:
        scattering_rate = [events/time for events, time in zip(scattering_data[scattering_type, :], time_spent)]
        all_scattering_rates.append(scattering_rate)
        # ax.plot (segments, data[scattering_type, :], '-')
        ax.plot (segments, scattering_rate, '-')
    ax.set_xlabel('Y (μm)', fontsize=12)
    ax.set_ylabel('Scattering rate (1/ns)', fontsize=12)
    legend = ["Sidewalls diffuse", "Sidewalls specular", "Top & bottom diffuse", "Top & bottom specular",
            "Holes diffuse", "Holes specular", "Hot side", "Internal", "Pillars diffuse", "Pillars specular"]
    ax.legend(legend, loc = 'upper right')
    plt.yscale('log')
    ax.set_ylim(bottom=1.0)
    fig.savefig("Scattering rates.pdf", dpi=300, format = 'pdf', bbox_inches="tight")
    if output_in_terminal: plt.show()

    # Save the file:
    header = "Y [um], Sidewalls diffuse [1/ns], Sidewalls specular [1/ns], Top & bottom diffuse [1/ns], Top & bottom specular [1/ns], Holes diffuse [1/ns], Holes specular [1/ns], Hot side [1/ns], Internal [1/ns], Pillars diffuse [1/ns], Pillars specular [1/ns]"
    # all_scattering_rates = np.transpose(all_scattering_rates)
    np.savetxt("Scattering rates.csv", np.vstack((segments, all_scattering_rates)).T, fmt='%1.2e', delimiter=",", header=header)

    # Ratio plot:
    # fig, ax = plt.subplots()
    # ax.plot (segments, data[0, :]/data[1, :], '-')
    # ax.set_xlabel('Y (μm)', fontsize=12)
    # ax.set_ylabel('Diffuse/scattering ratio', fontsize=12)
    # fig.savefig("Scattering statistics ratio.pdf", dpi=300, format = 'pdf', bbox_inches="tight")
    # if output_in_terminal: plt.show()


def output_time_spent_statistics(time_in_segments):
    '''Plot time that phonons spent in different segments'''
    # Calculate the segment coordinates:
    length_of_one_segment = length/number_of_length_segments
    segments = [(length_of_one_segment/2 + i*length_of_one_segment)*1e6 for i in range(number_of_length_segments)]

    # Calculate total time for all phonon groups and save it:
    total_time_in_segments = [np.sum(time_in_segments[segment, :])*1e9 for segment in range(number_of_length_segments)]
    np.savetxt('Time spent in segments.csv', np.vstack((segments, total_time_in_segments)).T, fmt='%1.3e', delimiter=",", header="Y [um], Time [ns]")

    # Plot the time spent in segments:
    fig, ax = plt.subplots()
    ax.plot (segments, total_time_in_segments, '-')
    ax.set_xlabel('Y (μm)', fontsize=12)
    ax.set_ylabel('Time spent (ns)', fontsize=12)
    fig.savefig("Time spent in segments.pdf", dpi=300, format = 'pdf', bbox_inches="tight")
    if output_in_terminal: plt.show()


def output_information(start_time, simulation_scheme):
    '''This function outputs the simulation information into the Information.txt file'''
    exit_angles = np.loadtxt("All exit angles.csv")
    percentage=int(100*np.count_nonzero(exit_angles)/(number_of_phonons+2*number_of_phonons*((simulation_scheme==2)*1)))
    print ("\n", percentage, '% of phonons reached the cold side')
    print ("The simulation took about", int((time.time()-start_time)//60), "min. to run")
    with open("Information.txt","w+") as f:
        info = (f'The simulation finished on {time.strftime("%d %B %Y")}, at {time.strftime("%H:%M")}.',
               f'\nIt took about {int((time.time()-start_time)//60)} min to run.\n',
               f'\nNumber of phonons = {number_of_phonons}',
               f'\nNumber of timesteps = {number_of_timesteps}',
               f'\nLength of a timestep = {timestep} s',
               f'\nTemperature = {T} K\n',
               f'\nLength = {length*1e9} nm',
               f'\nWidth = {width*1e9} nm',
               f'\nThickness = {thickness*1e9} nm\n',
               f'\nSide wall roughness = {side_wall_roughness*1e9} nm',
               f'\nHole roughness = {hole_roughness*1e9} nm',
               f'\nTop roughness = {top_roughness*1e9} nm',
               f'\nBottom roughness = {bottom_roughness*1e9} nm\n',
               f'\nLattice type = {hole_lattice_type}',
               f'\nPeriod in x direction = {period_x*1e9} nm',
               f'\nPeriod in y direction = {period_y*1e9} nm',
               f'\nDiameter of the holes = {circular_hole_diameter*1e9} nm',
               f'\nHorizontal dimension of the holes = {rectangular_hole_side_x*1e9} nm',
               f'\nVertical dimension of the holes = {rectangular_hole_side_y*1e9} nm\n',
               f'\n{percentage}% of phonons reached the cold side\n')
        f.writelines(info)

def output_general_statistics_on_scattering_events():
    '''This function calculated and outputs general statistics on scattering events'''
    stat = np.loadtxt("Statistics.csv", dtype='float')
    with open("Information.txt","a") as f:

        # Calculate the percentage of different scattering events:
        avg_scat = np.sum(stat)/number_of_phonons
        scat_on_walls = 100*(stat[0] + stat[1]) / np.sum(stat)
        scat_on_walls_diff = 100*stat[0]/(stat[0] + stat[1])
        scat_on_walls_spec = 100*stat[1]/(stat[0] + stat[1])
        scat_on_topbot = 100*(stat[2] + stat[3]) / np.sum(stat)
        scat_on_topbot_diff = 100*stat[2]/(stat[2] + stat[3])
        scat_on_topbot_spec = 100*stat[3]/(stat[2] + stat[3])
        if holes == 'yes':
            scat_on_holes = 100*(stat[4] + stat[5]) / np.sum(stat)
            scat_on_holes_diff = 100*stat[4]/(stat[4] + stat[5])
            scat_on_holes_spec = 100*stat[5]/(stat[4] + stat[5])
        if pillars == 'yes':
            scat_on_pillars = 100*(stat[8] + stat[9]) / np.sum(stat)
            scat_on_pillars_diff = 100*stat[8]/(stat[8] + stat[9])
            scat_on_pillars_spec = 100*stat[9]/(stat[8] + stat[9])
        retherm = 100*stat[6]/np.sum(stat)
        internal = 100*stat[7]/np.sum(stat)

        # Write it into the Information.txt
        f.writelines(['\nOn average, each phonon experienced %.2f scattering events' % avg_scat])
        f.writelines(['\n%.2f%% - scattering on side walls' % scat_on_walls,' (%.2f%% - diffuse,' % scat_on_walls_diff,' %.2f%% - specular)' % scat_on_walls_spec])
        f.writelines(['\n%.2f%% - scattering on top and bottom walls' % scat_on_topbot,' (%.2f%% - diffuse,' % scat_on_topbot_diff,' %.2f%% - specular)' % scat_on_topbot_spec])
        if holes == 'yes':
            f.writelines(['\n%.2f%% - scattering on hole walls' % scat_on_holes,' (%.2f%% - diffuse,' % scat_on_holes_diff,' %.2f%% - specular)' % scat_on_holes_spec])
        if pillars == 'yes':
            f.writelines(['\n%.2f%% - scattering on pillar walls' % scat_on_pillars,' (%.2f%% - diffuse,' % scat_on_pillars_diff,' %.2f%% - specular)' % scat_on_pillars_spec])
        f.writelines(['\n%.2f%% - rethermalization at the hot side' % retherm])
        f.writelines(['\n%.2f%% - internal scattering processes' % internal])
    os.remove("Statistics.csv")


def run_one_phonon(phonon_properties, statistics_of_scattering_events, maps_and_profiles, scattering_maps):
    '''This function runs one phonon through the system and returns various parametes of this run'''

    # Create some empty arrays and variables:
    x = np.zeros((number_of_timesteps))
    y = np.zeros((number_of_timesteps))
    z = np.zeros((number_of_timesteps))
    path_num = 0
    free_paths = [0.0]
    free_paths_along_y = [0.0]
    time_since_previous_scattering = 0.0
    exit_theta = 0.0
    travel_time = 0.0
    detected_frequency = 0.0

    # Get the initial properties and coordinates of the phonon at hot side:
    frequency, polarization, speed = phonon_properties[0:3]
    x[0], y[0], z[0], theta, phi = initialization()
    initial_theta = theta

    # If we use grey approximation, than expected time of internal scattering is simply mfp/speed:
    if use_gray_approximation_mfp:
        time_of_internal_scattering = gray_approximation_mfp/speed
    else:
        time_of_internal_scattering = internal_scattering_time_calculation(frequency, polarization)

    # Run the phonon step-by-step, until it reaches the hot side or the number of steps reaches maximum
    for i in range(1, number_of_timesteps):
        internal_scattering_type = 'none'
        reinitialization_scattering_type = 'none'

        # Check if phonon has not reached the cold side yet
        if phonon_is_in_system(x[i-1], y[i-1]):

            # Check if different scattering events happened during this time step:
            theta, phi, internal_scattering_type = internal_scattering(theta, phi, time_since_previous_scattering,
                    time_of_internal_scattering)

            theta, phi, surface_scattering_types = surface_scattering(x[i-1], y[i-1], z[i-1], theta, phi, frequency,
                    hole_coordinates, hole_shapes, pillar_coordinates,speed)

            theta, phi, reinitialization_scattering_type, x[i-1], y[i-1], z[i-1] = reinitialization(x[i-1],
                    y[i-1], z[i-1], theta, phi, speed)

            # If there was any scattering event, record it for future analysis:
            statistics_of_scattering_events = scattering_events_statistics_calculation(statistics_of_scattering_events,
                    surface_scattering_types, reinitialization_scattering_type, internal_scattering_type, y[i-1])

            # If there was no diffuse scattering event, then we keep measuring the paths and time:
            if (internal_scattering_type != 'diffuse') and (reinitialization_scattering_type != 'diffuse') \
                    and (all(i != 'diffuse' for i in surface_scattering_types)):

                # Increase the path and time since previous diffuse scattering by one timestep:
                free_paths[path_num] += speed*timestep
                time_since_previous_scattering += timestep

                # For serpentine lattice, we measure the free path along the structure in a special way:
                if hole_lattice_type == 'serpentine':
                    if abs(x[i-1]) < (width/2 - 155e-9):
                        free_paths_along_y[path_num] += speed*timestep*abs(cos(phi))*abs(sin(theta))
                    else:
                        free_paths_along_y[path_num] += speed*timestep*abs(cos(phi))*abs(cos(theta))
                # For all other lattices, we measure the free path along the structure along y:
                else:
                    free_paths_along_y[path_num] += speed*timestep*abs(cos(phi))*abs(cos(theta))

            # If diffuse scattering occurred, we reset phonon path and time since previous diffuse scattering:
            else:
                free_paths.append(0.0)
                free_paths_along_y.append(0.0)
                path_num += 1
                time_since_previous_scattering = 0.0
                # Calculate the expected time of internal scattering again (just for better randomization):
                if use_gray_approximation_mfp:
                    time_of_internal_scattering = gray_approximation_mfp/speed
                else:
                    time_of_internal_scattering = internal_scattering_time_calculation(frequency, polarization)

            # If we decided to record scattering map, here we update it, if there was any scattering:
            if output_scattering_map:
                scattering_maps = scattering_map_calculation(x[i-1], y[i-1], scattering_maps, internal_scattering_type,surface_scattering_types)
            maps_and_profiles = maps_and_profiles_calculation(x[i-1], y[i-1], maps_and_profiles, phonon_properties, i, theta, phi)

            # Phonon makes a step forward:
            x[i], y[i], z[i] = move(x[i-1], y[i-1], z[i-1], theta, phi, speed)

        # If the phonon reached cold side, record a few parameters and break the loop:
        else:
            exit_theta = theta          # Angle at the cold side
            travel_time = i*timestep    # Time it took to reach the cold side
            detected_frequency = frequency*(abs(x[i-1]) < frequency_detector_size/2.0)
            break

    flight_characteristics = [initial_theta, exit_theta, free_paths, free_paths_along_y, travel_time]
    return (flight_characteristics, x, y, z, statistics_of_scattering_events, maps_and_profiles,
            scattering_maps, detected_frequency)


def main1():
    '''This is the main function, which works under Debye approximation and should be used to simulate phonon paths at low temperatures'''
    print ('Simulation for', output_folder_name,' ')
    simulation_scheme = 1
    start_time = time.time()
    progress = -1

    # Create the folder if it does not exists and copy parameters.py there:
    if not os.path.exists(output_folder_name):
        os.makedirs(output_folder_name)
    shutil.copy('parameters.py', output_folder_name)

    # Create empty arrays
    all_initial_angles, all_exit_angles, all_free_paths, all_free_paths_along_y, all_frequencies, all_detected_frequencies, \
            all_group_velocities, all_travel_times = ([] for i in range(8))
    # statistics_of_scattering_events = [0]*10
    statistics_of_scattering_events = np.zeros((number_of_length_segments, 10))
    time_in_segments = np.zeros((number_of_length_segments, number_of_phonons//number_of_phonons_in_a_group))
    maps_and_profiles, scattering_maps = create_empty_maps()

    # For each group and for each phonon in this group:
    for i in range(number_of_phonons//number_of_phonons_in_a_group):    # To reduce RAM usage phonons are simulated in small groups
        x,y,z = (np.zeros((number_of_timesteps, number_of_phonons_in_a_group)) for index in range(3))
        for j in range(number_of_phonons_in_a_group):
            progress = progress_bar(i,j,progress,simulation_scheme)     # This is just to print the progress bar
            phonon_properties = phonon_properties_assignment()          # Get initial phonon properties: frequency, polarization, and speed
            phonon_properties.append(i*number_of_phonons_in_a_group+j)  # Also add phonon number to phonon properties

            # Run this phonon through the structure:
            flight_characteristics, x[:,j], y[:,j], z[:,j], statistics_of_scattering_events, maps_and_profiles, \
                    scattering_maps,detected_frequency = run_one_phonon(phonon_properties,
                    statistics_of_scattering_events, maps_and_profiles, scattering_maps)

            # Record the properties returned for this phonon:
            all_initial_angles.append(flight_characteristics[0])
            all_exit_angles.append(flight_characteristics[1])
            all_free_paths.extend(flight_characteristics[2])
            all_free_paths_along_y.extend(flight_characteristics[3])
            all_travel_times.append(flight_characteristics[4])
            all_frequencies.append(phonon_properties[0])
            all_detected_frequencies.append(detected_frequency)
            all_group_velocities.append(phonon_properties[2])
        time_in_segments = record_time_in_segments(time_in_segments, y, i)

    # Write files and make various plots:
    write_files(all_free_paths, all_free_paths_along_y, all_frequencies, all_exit_angles, all_initial_angles,
            all_group_velocities, statistics_of_scattering_events, all_travel_times, all_detected_frequencies)
    output_distributions()
    output_thermal_map(maps_and_profiles[0])
    if output_scattering_map:
        output_scattering_maps(scattering_maps)
    output_profiles(maps_and_profiles)
    output_thermal_conductivity(maps_and_profiles)
    output_trajectories(x,y,z,number_of_phonons_in_a_group)
    output_information(start_time, simulation_scheme)
    output_general_statistics_on_scattering_events()
    output_time_spent_statistics(time_in_segments)
    output_scattering_statistics(statistics_of_scattering_events)
    #coordinates=np.zeros((number_of_timesteps,number_of_phonons_in_a_group*3))


def main2():
    '''This is the main function, which calculates thermal conductivity by integrating bulk dispersion'''
    print ('Simulation for',output_folder_name,'has started')
    simulation_scheme = 2
    start_time = time.time()
    progress = -1
    all_initial_angles,all_exit_angles,all_free_paths,all_free_paths_along_y,all_frequencies,all_group_velocities,all_travel_times = ([] for i in range(8))
    statistics_of_scattering_events = [0]*10
    maps_and_profiles, scattering_maps = create_empty_maps()

    # Create the folder if it does not exists and copy parameters.py there:
    if not os.path.exists(output_folder_name):
        os.makedirs(output_folder_name)
    shutil.copy('parameters.py', output_folder_name)

    mean_free_path = np.zeros((number_of_phonons))
    thermal_conductivity = 0
    cummulative_conductivity = np.zeros((number_of_phonons,6))
    for branch in range(3):
        x,y,z = (np.zeros((number_of_timesteps,number_of_phonons_in_a_group)) for i in range(3))
        for j in range(0,number_of_phonons):
            progress = progress_bar(branch,j,progress,simulation_scheme)
            phonon_properties, w, K, dK  = phonon_properties_assignment_2(j,branch)
            flight_characteristics,x[:,j],y[:,j],z[:,j],statistics_of_scattering_events,maps_and_profiles, scattering_maps = run_one_phonon(phonon_properties,statistics_of_scattering_events,maps_and_profiles,scattering_maps)

            # Record the properties returned for this phonon:
            all_initial_angles.append(flight_characteristics[0])
            all_exit_angles.append(flight_characteristics[1])
            all_free_paths.extend(flight_characteristics[2])
            all_free_paths_along_y.extend(flight_characteristics[3])
            all_travel_times.append(flight_characteristics[4])
            all_frequencies.append(phonon_properties[0])
            all_group_velocities.append(phonon_properties[2])

            # Calculate the thermal conductivity:
            frequency, polarization, speed = phonon_properties
            mean_free_path[j]=sum(flight_characteristics[2])/len(flight_characteristics[2])                         # Average of all free paths for this phonon
            heat_capacity=k*((hbar*w/(k*T))**2)*exp(hbar*w/(k*T))/((exp(hbar*w/(k*T))-1)**2)                        # Ref. PRB 88 155318 (2013)
            thermal_conductivity+=(1/(6*(pi**2)))*heat_capacity*(speed**2)*(mean_free_path[j]/speed)*(K**2)*dK      # Eq.3 from Physical Review 132 2461 (1963)

            cummulative_conductivity[j,branch]=speed/frequency
            cummulative_conductivity[j,branch+3]=(1/(6*(pi**2)))*heat_capacity*(speed**2)*(mean_free_path[j]/speed)*(K**2)*dK

    # Write files and make various plots:
    write_files(all_free_paths,all_free_paths_along_y,all_frequencies,all_exit_angles,all_initial_angles,all_group_velocities,statistics_of_scattering_events,all_travel_times)
    output_trajectories(x,y,z,number_of_phonons)
    output_thermal_map(thermal_map)
    output_distributions()
    output_information(start_time, simulation_scheme)
    output_general_statistics_on_scattering_events()
    print ('Thermal conductivity =', thermal_conductivity)

    plt.figure(12)
    for i in range(3):
        plt.loglog (cummulative_conductivity[:,i]*1e9,cummulative_conductivity[:,i+3])
    plt.show()
    np.savetxt('Distribution of wavelengths.csv', cummulative_conductivity, fmt='%1.3e', delimiter=",")


if __name__ == "__main__":

    # Mode 1 is the main mode.
    if simulation_mode == 1:
        main1()

    # Mode 2 is basically Callaway-Holland model but the relaxation time is actually measured
elif simulation_mode == 2:
    main2()
