import numpy as np
from random import random, choice
from math import pi, cos, sin, tan, exp, log, sqrt, atan, asin, acos
from scipy.constants import k, hbar

from parameters import *

def initialization():
    '''Assign initial position and angles of new phonons at the hot side'''

    # By default, phonons start at random position at the hot side:
    x = 0.49*width*(2*random()-1)
    y = 1e-12                       # It is practically zero
    z = 0.49*thickness*(2*random()-1)

    # In some lattices, phonons start from specific places:
    if hole_lattice_type == 'serpentine':
        x = -width/2+(155e-9)/2 + 0.4*(155e-9)*(2*random()-1)
    if hole_lattice_type == 'diode_with_wires':
        x = 0.4*(period_x-rectangular_hole_side_x)*(2*random()-1)
    if hole_lattice_type == 'turn':
        x = 0.4*(period_x*5)*(2*random()-1)
    if hole_lattice_type == 'turn90':
        x = 0.4*width/2.0*(2*random()-1) - period_x*3.5
    if hole_lattice_type == 'directional_source':
        x = (2*random()-1)*50e-9
    if hole_lattice_type == 'tesla valve forward':
        x = (2*random()-1)*10e-9 - 100e-9
    if hole_lattice_type == 'tesla valve backward':
        x = (2*random()-1)*10e-9 - 100e-9
    if hole_lattice_type == 'diamond_particle':
        x = 0.40*thickness*(2*random()-1)


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

    if material == 'Si':                                        # Ref. APL 95 161901 (2009)
        dispersion[:,0] = [k*12e9/(N-1) for k in range(N)]                                                              # Wavevectors
        dispersion[:,1] = [abs(1369.42*k-2.405e-8*(k**2)-9.70e-19*(k**3)) for k in dispersion[:,0]]                     # LA branch
        dispersion[:,2] = [abs(1081.74*k-7.711e-8*(k**2)+5.674e-19*(k**3)+7.967e-29*(k**4)) for k in dispersion[:,0]]   # TA branch
        dispersion[:,3] = dispersion[:,2]                                                                               # TA branch

    if material == 'SiC':                                       # https://journals.aps.org/prb/pdf/10.1103/PhysRevB.50.17054
        dispersion[:,0] = [k*14414281503/(N-1) for k in range(N)]                                                       # Wavevectors
        dispersion[:,1] = [abs(-3.48834e-18*(k**3)+1.7604452e-08*(k**2)+1737.36296*k) for k in dispersion[:,0]]         # LA branch
        dispersion[:,2] = [abs(-2.21696e-19*(k**3)-3.4366886e-08*(k**2)+1077.98941*k) for k in dispersion[:,0]]         # TA branch
        dispersion[:,3] = dispersion[:,2]                                                                               # TA branch

    if material == 'Diamond':                                   # https://www.sciencedirect.com/science/article/pii/S0008622315003358
        dispersion[:,0] = [k*11707071561.7/(N-1) for k in range(N)]                                                     # Wavevectors
        dispersion[:,1] = [abs(-1.347265e-18*(k**3)-8.855338e-08*(k**2)+4309.95222*k) for k in dispersion[:,0]]         # LA branch
        dispersion[:,2] = [abs(-5.042335e-18*(k**3)-4.104260e-08*(k**2)+3185.66561*k) for k in dispersion[:,0]]         # TA branch
        dispersion[:,3] = dispersion[:,2]                                                                               # TA branch

    return dispersion


def phonon_properties_assignment():
    '''Assign phonon frequency (f) according to the Plank distribution at a given temperature T,
    chose polarization, and calculate group velocity from bulk dispersion'''
    if material == 'Si':
        default_speed = 6000                                          # [m/s] This is the speed for Debye approximation
    if material == 'SiC':
        default_speed = 6500                                          # [m/s] Need to change this probably!
    if material == 'Diamond':
        default_speed = 20000                                          # [m/s]
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


def create_empty_maps():
    '''This function just creates empty maps'''
    # Create empty profile map arrays:
    thermal_map = np.zeros((number_of_pixels_y,number_of_pixels_x))
    heat_flux_profile_x = np.zeros((number_of_pixels_x,number_of_timeframes))
    heat_flux_profile_y = np.zeros((number_of_pixels_y,number_of_timeframes))
    temperature_profile_x = np.zeros((number_of_pixels_x,number_of_timeframes))
    temperature_profile_y = np.zeros((number_of_pixels_y,number_of_timeframes))

    # Create empty scattering maps:
    diffuse_scattering_map_x  = []
    diffuse_scattering_map_y  = []
    specular_scattering_map_x = []
    specular_scattering_map_y = []
    internal_scattering_map_x = []
    internal_scattering_map_y = []

    # Package everything into one tuple:
    maps_and_profiles = [thermal_map,heat_flux_profile_x,heat_flux_profile_y,temperature_profile_x,temperature_profile_y]
    scattering_maps = [diffuse_scattering_map_x, diffuse_scattering_map_y, specular_scattering_map_x,
            specular_scattering_map_y, internal_scattering_map_x, internal_scattering_map_y]
    return maps_and_profiles, scattering_maps
