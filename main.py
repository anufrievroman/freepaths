# General libraries:
import os
import time
import sys
import shutil
import matplotlib.pyplot as plt
import numpy as np

# Specific functions:
from math import pi, cos, sin, exp
from scipy.constants import k, hbar

# Modules:
from parameters import *
from initialization import *
from output import *
from scattering import *
from lattices import *
from move import *

plt.style.use('plotstyle.mplstyle')


def progress_bar(i, j, old_progress, scheme):
    '''This is a progress bar that outputs the progress each one percent but not more often'''
    if scheme == 1:
        progress = 100*(i*number_of_phonons_in_a_group+j)//number_of_phonons
    elif scheme == 2:
        progress = 100*(i*number_of_phonons+j)//(number_of_phonons*3)

    # Only update the progress bar if it changed by 1%
    if progress > old_progress:
        sys.stdout.write('\r' + 'Progress: ' + str(progress) + '%')
        sys.stdout.flush()
    return progress


def phonon_is_in_system(x, y):
    '''This function checks if the phonon at this timestep is still in the system and did not reach the cold side.
    Depending on where we set the cold side, we check if phonon crossed that line'''

    small_offset = 10e-9
    if cold_side == 'top':
        phonon_is_in_system = (y < length and y > 0)
    elif cold_side == 'right':
        phonon_is_in_system = ((y < length-1.1e-6) or (y > length-1.1e-6 and x < width/2.0-small_offset))
    elif cold_side == 'top and right':
        phonon_is_in_system = ((y < 1.0e-6) or (y > 1.0e-6 and x < width/2.0-small_offset and y < length))
    return phonon_is_in_system


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

    # Run the phonon step-by-step, until it reaches the hot side or the number of steps reaches maximum:
    for i in range(1, number_of_timesteps):
        internal_scattering_type = 'none'
        reinitialization_scattering_type = 'none'

        # Check if phonon has not reached the cold side yet:
        if phonon_is_in_system(x[i-1], y[i-1]):

            # Check if different scattering events happened during this time step:
            theta, phi, internal_scattering_type = internal_scattering(theta, phi, time_since_previous_scattering,
                    time_of_internal_scattering)

            theta, phi, surface_scattering_types = surface_scattering(x[i-1], y[i-1], z[i-1], theta, phi, frequency,
                    hole_coordinates, hole_shapes, pillar_coordinates,speed)

            # theta, phi, reinitialization_scattering_type, x[i-1], y[i-1], z[i-1] = reinitialization(x[i-1],
                    # y[i-1], z[i-1], theta, phi, speed)

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
                scattering_maps = scattering_map_calculation(x[i-1], y[i-1], scattering_maps,
                                                            internal_scattering_type,surface_scattering_types)
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
        x, y, z = (np.zeros((number_of_timesteps, number_of_phonons_in_a_group)) for index in range(3))
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
    output_trajectories(x, y, z, number_of_phonons_in_a_group)
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
    all_initial_angles, all_exit_angles, all_free_paths, all_free_paths_along_y, all_frequencies, all_detected_frequencies, \
            all_group_velocities, all_travel_times = ([] for i in range(8))
    statistics_of_scattering_events = np.zeros((number_of_length_segments, 10))
    maps_and_profiles, scattering_maps = create_empty_maps()

    # Create the folder if it does not exists and copy parameters.py there:
    if not os.path.exists(output_folder_name):
        os.makedirs(output_folder_name)
    shutil.copy('parameters.py', output_folder_name)

    mean_free_path = np.zeros((number_of_phonons))
    thermal_conductivity = 0
    cummulative_conductivity = np.zeros((number_of_phonons,6))
    for branch in range(3):
        x, y, z = (np.zeros((number_of_timesteps, number_of_phonons_in_a_group)) for i in range(3))
        for j in range(0,number_of_phonons):
            progress = progress_bar(branch,j,progress,simulation_scheme)
            phonon_properties, w, K, dK  = phonon_properties_assignment_2(j,branch)

            # Run this phonon through the structure:
            flight_characteristics, x[:,j], y[:,j], z[:,j], statistics_of_scattering_events, maps_and_profiles, \
                    scattering_maps, detected_frequency = run_one_phonon(phonon_properties,
                    statistics_of_scattering_events, maps_and_profiles, scattering_maps)

            # Record the properties returned for this phonon:
            all_initial_angles.append(flight_characteristics[0])
            all_exit_angles.append(flight_characteristics[1])
            all_free_paths.extend(flight_characteristics[2])
            all_free_paths_along_y.extend(flight_characteristics[3])
            all_travel_times.append(flight_characteristics[4])
            all_frequencies.append(phonon_properties[0])
            all_group_velocities.append(phonon_properties[2])
            all_detected_frequencies.append(detected_frequency)

            # Calculate the thermal conductivity:
            frequency, polarization, speed = phonon_properties
            mean_free_path[j]=sum(flight_characteristics[2])/len(flight_characteristics[2])  # Average of all free paths for this phonon
            heat_capacity=k*((hbar*w/(k*T))**2)*exp(hbar*w/(k*T))/((exp(hbar*w/(k*T))-1)**2) # Ref. PRB 88 155318 (2013)
            thermal_conductivity+=(1/(6*(pi**2)))*heat_capacity*(speed**2)*(mean_free_path[j]/speed)*(K**2)*dK  # Eq.3 from Physical Review 132 2461 (1963)

            cummulative_conductivity[j,branch]=speed/frequency
            cummulative_conductivity[j,branch+3]=(1/(6*(pi**2)))*heat_capacity*(speed**2)*(mean_free_path[j]/speed)*(K**2)*dK

    # Write files and make various plots:
    write_files(all_free_paths, all_free_paths_along_y, all_frequencies, all_exit_angles, all_initial_angles,
            all_group_velocities, statistics_of_scattering_events, all_travel_times, all_detected_frequencies)
    output_trajectories(x, y, z, number_of_phonons)
    output_distributions()
    output_thermal_map(maps_and_profiles[0])
    output_information(start_time, simulation_scheme)
    output_general_statistics_on_scattering_events()
    print ('Thermal conductivity =', thermal_conductivity)

    plt.figure(12)
    for i in range(3):
        plt.loglog (cummulative_conductivity[:,i]*1e9,cummulative_conductivity[:,i+3])
    # plt.show()
    np.savetxt('Distribution of wavelengths.csv', cummulative_conductivity, fmt='%1.3e', delimiter=",")


if __name__ == "__main__":

    # Mode 1 is the main mode.
    if simulation_mode == 1:
        main1()

    # Mode 2 is basically Callaway-Holland model but the relaxation time is actually measured
    elif simulation_mode == 2:
        main2()
