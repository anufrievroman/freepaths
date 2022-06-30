import time
import random
import sys
import os
from math import pi, cos
from scipy.constants import hbar

import numpy as np
import matplotlib.pyplot as plt

from parameters import *


def distribution_calculation(filename, data_range, number_of_nodes):
    '''This function calculates distribution of numbers (histogram) in a given file'''
    data = np.loadtxt(filename)
    if data_range is None:
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

    # Scattering on pillars:
    statistics_of_scattering_events[segment, 8] += 1 if surface_scattering_types[3] == 'diffuse' else 0
    statistics_of_scattering_events[segment, 9] += 1 if surface_scattering_types[3] == 'specular' else 0
    return statistics_of_scattering_events


def scattering_map_calculation(x,y,scattering_maps,internal_scattering_type,surface_scattering_types):
    '''This function records the place where a scattering event occurred according to the event type'''
    # Diffuse surface scattering:
    if any(scattering == 'diffuse' for scattering in surface_scattering_types):
        scattering_maps[0].append(x)
        scattering_maps[1].append(y)

    # Specular surface scattering:
    elif any(scattering == 'specular' for scattering in surface_scattering_types):
        scattering_maps[2].append(x)
        scattering_maps[3].append(y)

    # Internal scattering:
    elif internal_scattering_type == 'diffuse':
        scattering_maps[4].append(x)
        scattering_maps[5].append(y)
    return scattering_maps


def record_time_in_segments(time_in_segments, y, group):
    '''Record how long phonons stay in different segments'''
    for segment_number in range(number_of_length_segments):
        segment_begining = segment_number*(length/number_of_length_segments)
        segment_end = (segment_number+1)*(length/number_of_length_segments)
        for phonon in range(number_of_phonons_in_a_group):
            times = [timestep for point in y[:, phonon] if ((segment_begining <= point < segment_end) and point != 0.0)]
            time_in_segments[segment_number, group] += sum(times)
    return time_in_segments


def maps_and_profiles_calculation(x,y,maps_and_profiles,phonon_properties,timestep_number,theta,phi):
    '''This function registers the phonon in the pixel corresponding to its current position
    and at certain timesteps and calculates thermal maps and thermal profiles along different axes'''
    # Unpacking stuff from the tuple container:
    thermal_map, heat_flux_profile_x, heat_flux_profile_y, temperature_profile_x, temperature_profile_y = maps_and_profiles
    frequency, _, speed = phonon_properties[0:3]

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
    elif material == 'Diamond':
        material_density = 3500 # [kg/m^3]

    # Here we arbitraraly correct the volume of the unit cells:
    if pillars=='yes':
        Vcell_x+=2.5*0.3333*pillar_height*(circular_hole_diameter/2)**2
        Vcell_y+=2.5*0.3333*pillar_height*(circular_hole_diameter/2)**2

    # Prevent error if the phonon is outside the structure:
    if (0 <= index_x < number_of_pixels_x) and (0 <= index_y < number_of_pixels_y):
        # Record energy h*w of this phonon into the pixel of thermal map:
        thermal_map[index_y,index_x] += hbar*2*pi*frequency

        # Record energy of this phonon into flux and temperature profiles:
        timeframe_number = int(((timestep_number+random.randint(0, number_of_timesteps))*timestep*number_of_timeframes) // (number_of_timesteps*timestep))
        if timeframe_number < number_of_timeframes:
            heat_flux_profile_x[index_x,timeframe_number]   += hbar*2*pi*frequency*cos(theta)*abs(cos(phi))*speed/Vcell_x
            heat_flux_profile_y[index_y,timeframe_number]   += hbar*2*pi*frequency*cos(theta)*abs(cos(phi))*speed/Vcell_y
            temperature_profile_x[index_x,timeframe_number] += hbar*2*pi*frequency/(specific_heat_capacity*material_density)/Vcell_x
            temperature_profile_y[index_y,timeframe_number] += hbar*2*pi*frequency/(specific_heat_capacity*material_density)/Vcell_y

    # Pack everything back into the tuple container and return:
    maps_and_profiles = [thermal_map, heat_flux_profile_x, heat_flux_profile_y, temperature_profile_x, temperature_profile_y]
    return maps_and_profiles


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
            data[j, i] = value                            # There is some error, where j comes from?
    np.savetxt("Scattering map.csv", data, fmt='%1.2e', delimiter=",")


def output_thermal_conductivity(maps_and_profiles):
    '''Calculate the thermal conductivity for each time interval from heat flux
    and temperature profiles accumulated in that interval'''
    _, J_profiles_x, J_profiles_y, T_profiles_x, T_profiles_y = maps_and_profiles
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
