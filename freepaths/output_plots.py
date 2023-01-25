"""Module that calculates and outputs vaious plots and distributions from the saved files"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from freepaths.config import cf
import matplotlib.pyplot as plt

# Style of the plots:
plt.rcParams['font.family'] = "Arial"
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['lines.linewidth'] = 1.0
plt.rcParams['xtick.direction'] = "in"
plt.rcParams['ytick.direction'] = "in"
plt.rcParams['xtick.minor.visible'] = False
plt.rcParams['ytick.minor.visible'] = False
plt.rcParams['legend.frameon'] = False
plt.rcParams['figure.autolayout'] = True
plt.rcParams['figure.figsize'] = [5, 3.5]
plt.rcParams['figure.dpi'] = 200
plt.rcParams['savefig.dpi'] = 200
plt.rcParams['legend.fontsize'] = 8


def distribution_calculation(filename, data_range, number_of_nodes):
    """Calculate distribution of numbers (histogram) in a given file"""
    data = np.loadtxt(filename, encoding='utf-8')
    if data_range is None:
        data_range = np.max(data)
    distribution = np.zeros((number_of_nodes, 2))
    distribution[:, 0] = np.linspace(0, data_range, number_of_nodes)
    distribution[:, 1], _ = np.histogram(data[data != 0], number_of_nodes, range=(0, data_range))
    return distribution


def angle_distribution_calculation():
    """Analyse measured phonon angles and create their distribution"""
    all_exit_angles = np.loadtxt("Data/All exit angles.csv", dtype='float', encoding='utf-8')
    initial_angles = np.loadtxt("Data/All initial angles.csv", dtype='float', encoding='utf-8')
    distribution = np.zeros((180, 3))
    distribution[:, 0] = range(-90, 90)
    exit_angles = all_exit_angles[all_exit_angles != 0]
    distribution[:, 1], _ = np.histogram(np.degrees(exit_angles), 180, range=(-90, 90))
    distribution[:, 2], _ = np.histogram(np.degrees(initial_angles), 180, range=(-90, 90))
    return distribution


def wavelength_distribution_calculation(number_of_nodes):
    """Calculate phonon wavelength distribution from their frequencies and velocities"""
    frequencies = np.loadtxt("Data/All initial frequencies.csv", encoding='utf-8')
    speeds = np.loadtxt("Data/All group velocities.csv", encoding='utf-8')
    wavelengths = np.zeros((len(speeds)))
    wavelengths[:] = speeds[:] / frequencies[:]
    data_range = np.amax(wavelengths)
    distribution = np.zeros((number_of_nodes, 2))
    distribution[:, 0] = np.linspace(0, data_range, number_of_nodes)
    distribution[:, 1], _ = np.histogram(wavelengths, number_of_nodes, range=(0, data_range))
    return distribution


def plot_angle_distribution():
    """Plot distribution of angles"""
    angle_distributions = angle_distribution_calculation()
    fig, ax = plt.subplots()
    ax.plot(angle_distributions[:, 0], angle_distributions[:, 1], 'b')
    ax.plot(angle_distributions[:, 0], angle_distributions[:, 2], 'r')
    ax.set_xlabel('Angle (degree)', fontsize=12)
    ax.set_ylabel('Number of phonons', fontsize=12)
    ax.legend(["At hot side", "At cold side"])
    fig.savefig("Distribution of angles.pdf", dpi=300, format='pdf', bbox_inches="tight")
    if cf.plots_in_terminal: plt.show()
    np.savetxt('Data/Distribution of angles.csv', angle_distributions, fmt='%1.3e', delimiter=",")


def plot_free_path_distribution():
    """Plot distribution of free path"""
    filename = "Data/All free paths.csv"
    free_path_distribution = distribution_calculation(filename, None, cf.number_of_nodes)
    fig, ax = plt.subplots()
    ax.plot(free_path_distribution[:, 0] * 1e6, free_path_distribution[:, 1])
    ax.set_xlabel('Free flights (μm)', fontsize=12)
    ax.set_ylabel('Number of flights', fontsize=12)
    # ax.set_xlim([0, max(free_path_distribution[:,0])*1e6])
    fig.savefig("Distribution of free paths.pdf", dpi=300, format='pdf', bbox_inches="tight")
    if cf.plots_in_terminal: plt.show()
    np.savetxt('Data/Distribution of free paths.csv', free_path_distribution, fmt='%1.3e', delimiter=",")


def plot_frequency_distribution():
    """Plot distribution of frequencies"""
    filename = "Data/All initial frequencies.csv"
    frequency_distribution = distribution_calculation(filename, None, cf.number_of_nodes)
    fig, ax = plt.subplots()
    ax.plot(frequency_distribution[:, 0], frequency_distribution[:, 1])
    ax.set_xlabel('Frequency (Hz)', fontsize=12)
    ax.set_ylabel('Number of phonons', fontsize=12)
    fig.savefig("Distribution of initial frequencies.pdf", dpi=300, format='pdf', bbox_inches="tight")
    if cf.plots_in_terminal: plt.show()
    np.savetxt('Data/Distribution of initial frequencies.csv', frequency_distribution, fmt='%1.3e', delimiter=",")


def plot_wavelength_distribution():
    """Plot distribution of wavelength"""
    wavelength_distribution = wavelength_distribution_calculation(cf.number_of_nodes)
    fig, ax = plt.subplots()
    ax.plot(wavelength_distribution[:, 0] * 1e9, wavelength_distribution[:, 1])
    ax.set_xlabel('Wavelength (nm)', fontsize=12)
    ax.set_ylabel('Number of phonons', fontsize=12)
    fig.savefig("Distribution of wavelengths.pdf", dpi=300, format='pdf', bbox_inches="tight")
    if cf.plots_in_terminal: plt.show()
    np.savetxt('Data/Distribution of wavelengths.csv', wavelength_distribution, fmt='%1.3e', delimiter=",")


def plot_travel_time_distribution():
    """Plot distribution of wavelength"""
    travel_time_distribution = distribution_calculation("Data/All travel times.csv", None, cf.number_of_nodes)
    fig, ax = plt.subplots()
    ax.plot(travel_time_distribution[:, 0] * 1e9, travel_time_distribution[:, 1])
    ax.set_xlabel('Travel time (ns)', fontsize=12)
    ax.set_ylabel('Number of phonons', fontsize=12)
    fig.savefig("Distribution of travel times.pdf", dpi=300, format='pdf', bbox_inches="tight")
    if cf.plots_in_terminal: plt.show()
    np.savetxt('Data/Distribution of travel times.csv', travel_time_distribution, fmt='%1.3e', delimiter=",")


def plot_detected_frequency_distribution():
    """Plot distribution of detected frequencies"""
    detected_frequency_distribution = distribution_calculation("Data/All detected frequencies.csv", None, cf.number_of_nodes)
    fig, ax = plt.subplots()
    ax.plot(detected_frequency_distribution[:, 0], detected_frequency_distribution[:, 1])
    ax.set_xlabel('Frequency (Hz)', fontsize=12)
    ax.set_ylabel('Number of phonons', fontsize=12)
    fig.savefig("Distribution of detected frequencies.pdf", dpi=300, format='pdf', bbox_inches="tight")
    if cf.plots_in_terminal: plt.show()
    np.savetxt('Data/Distribution of detected frequencies.csv', detected_frequency_distribution, fmt='%1.3e', delimiter=",")


def plot_velocity_distribution():
    """Plot distribution of group velocities"""
    fig, ax = plt.subplots()
    speeds = np.loadtxt("Data/All group velocities.csv")
    frequencies = np.loadtxt("Data/All initial frequencies.csv")
    ax.plot(frequencies, speeds, '.')
    ax.set_xlabel('Frequency (Hz)', fontsize=12)
    ax.set_ylabel('Group velocity (m/s)', fontsize=12)
    fig.savefig('Group velocities.pdf', dpi=300, format='pdf', bbox_inches="tight")
    if cf.plots_in_terminal: plt.show()


def plot_time_in_segments():
    """Plot time spent in segments"""
    fig, ax = plt.subplots()
    segment, time = np.genfromtxt("Data/Time spent in segments.csv", unpack=True, delimiter=',', usecols=(0, 1), skip_header=1)
    ax.plot(segment, time, '-')
    ax.set_xlabel('Y (μm)', fontsize=12)
    ax.set_ylabel('Time spent (ns)', fontsize=12)
    fig.savefig("Time spent in segments.pdf", dpi=300, format='pdf', bbox_inches="tight")
    if cf.plots_in_terminal: plt.show()


def plot_thermal_conductivity():
    """Plot thermal conductivity against time segment"""
    fig, ax = plt.subplots()
    time, thermal_conductivity = np.genfromtxt("Data/Thermal conductivity.csv", unpack=True, delimiter=',', usecols=(0, 1), skip_header=1)
    ax.plot(time, thermal_conductivity, linewidth=1)
    ax.set_ylabel('Thermal conductivity (W/mK)', fontsize=12)
    ax.set_xlabel('Time (ns)', fontsize=12)
    fig.savefig("Thermal conductivity.pdf", dpi=300, format='pdf', bbox_inches="tight")
    if cf.plots_in_terminal: plt.show()


def plot_temperature_profile():
    """Plot profile of temperature for each time segment"""
    fig, ax = plt.subplots()
    data = np.genfromtxt("Data/Temperature profiles y.csv", unpack=True, delimiter=',', skip_header=1, encoding='utf-8')
    for timeframe in range(len(data) - 1):
        ax.plot(data[0], data[timeframe + 1], linewidth=1)
    ax.set_xlabel('Y (μm)', fontsize=12)
    ax.set_ylabel('Temperature (K)', fontsize=12)
    fig.savefig("Temperature profile.pdf", dpi=300, format='pdf', bbox_inches="tight")
    if cf.plots_in_terminal: plt.show()


def plot_heat_flux_profile():
    """Plot profile of heat flux for each time segment"""
    fig, ax = plt.subplots()
    data = np.genfromtxt("Data/Heat flux profiles y.csv", unpack=True, delimiter=',', skip_header=1, encoding='utf-8')
    for timeframe in range(len(data) - 1):
        ax.plot(data[0][1:], data[timeframe + 1][1:], linewidth=1)
    ax.set_xlabel('Y (μm)', fontsize=12)
    ax.set_ylabel('Heat flux (W/m^2)', fontsize=12)
    fig.savefig("Heat flux profile.pdf", dpi=300, format='pdf', bbox_inches="tight")
    if cf.plots_in_terminal: plt.show()


def plot_thermal_map():
    """Plot thermal map as color map"""
    fig = plt.figure()
    thermal_map = np.genfromtxt("Data/Thermal map.csv", unpack=False, delimiter=',', skip_header=0, encoding='utf-8')
    thermal_map = np.flipud(thermal_map)
    minimum_of_colorbar = 1e-20  # Cannot be zero!
    boundaries = [(-cf.width / 2) * 1e6, (cf.width / 2) * 1e6, 0, cf.length * 1e6]
    plt.imshow(thermal_map, cmap='hot', interpolation='none', extent=boundaries,
               norm=LogNorm(vmin=minimum_of_colorbar, vmax=np.amax(thermal_map)))
    plt.xlabel('X (μm)', fontsize=12)
    plt.ylabel('Y (μm)', fontsize=12)
    cbar = plt.colorbar()
    cbar.set_label('Energy density', rotation=90)
    fig.savefig("Thermal map.pdf", bbox_inches="tight")
    if cf.plots_in_terminal: plt.show()


def plot_scattering_map():
    """Plot the map of scattering events"""
    fig, ax = plt.subplots()
    filename = "Data/Scattering map.csv"
    spec_x, spec_y, diff_x, diff_y, int_x, int_y = np.genfromtxt(filename, unpack=True, delimiter=',', skip_header=1, encoding='utf-8')
    ax.plot(diff_x[diff_x != 0], diff_y[diff_y != 0], 'o', color='b', markersize=0.1, alpha=0.3)
    ax.plot(spec_x[spec_x != 0], spec_y[spec_y != 0], 'o', color='g', markersize=0.1, alpha=0.3)
    ax.plot(int_x[int_x != 0], int_y[int_y != 0], 'o', color='r', markersize=0.1, alpha=0.3)
    ax.set_xlabel('X (μm)', fontsize=12)
    ax.set_ylabel('Y (μm)', fontsize=12)
    ax.legend(["Diffuse", "Specular", "Internal"])
    ax.set_aspect('equal', 'datalim')
    fig.savefig("Scattering map.pdf", bbox_inches="tight")
    if cf.plots_in_terminal: plt.show()


def plot_trajectories():
    """Plot the phonon trajectories"""

    data = np.genfromtxt("Data/Phonon paths.csv", unpack=False, delimiter=',', skip_header=1, encoding='utf-8')

    # Create XY plot:
    fig, ax = plt.subplots()
    for index in range(cf.output_trajectories_of_first):
        x_coordinates = np.trim_zeros(data[:, 3 * index], trim='b')
        y_coordinates = np.trim_zeros(data[:, 3 * index + 1], trim='b')
        ax.plot(x_coordinates, y_coordinates, linewidth=0.2)
    ax.set_xlabel('X (μm)', fontsize=12)
    ax.set_ylabel('Y (μm)', fontsize=12)
    ax.set_aspect('equal', 'datalim')
    fig.savefig("Phonon paths XY.pdf", dpi=600, format='pdf', bbox_inches="tight")
    if cf.plots_in_terminal: plt.show()

    # Create YZ plot:
    fig, ax = plt.subplots()
    for index in range(cf.output_trajectories_of_first):
        y_coordinates = np.trim_zeros(data[:, 3 * index + 1], trim='b')
        z_coordinates = np.trim_zeros(data[:, 3 * index + 2], trim='b')
        ax.plot(y_coordinates, z_coordinates, linewidth=0.2)
    ax.set_xlabel('Y (μm)', fontsize=12)
    ax.set_ylabel('Z (μm)', fontsize=12)
    ax.set_aspect('equal', 'datalim')
    fig.savefig("Phonon paths YZ.pdf", dpi=300, format='pdf', bbox_inches="tight")
    if cf.plots_in_terminal: plt.show()


def plot_scattering_statistics():
    """Calculate and plot rates of different scattering events in each length segment"""
    # Load the data from files:
    filename = "Data/Scattering events statistics.csv"
    scattering_data = np.genfromtxt(filename, unpack=True, delimiter=',', skip_header=1, encoding='utf-8')
    filename = "Data/Time spent in segments.csv"
    segments, time_spent = np.genfromtxt(filename, unpack=True, delimiter=',', skip_header=1, encoding='utf-8')

    # Create the plot:
    fig, ax = plt.subplots()
    all_scattering_rates = []
    for scattering_type in range(scattering_data.shape[0]):
        # Calculate and plot scattering events per second:
        scattering_rate = [events / time for events, time in zip(scattering_data[scattering_type, :], time_spent)]
        all_scattering_rates.append(scattering_rate)
        ax.plot(segments, scattering_rate, '-')

    ax.set_xlabel('Y (μm)', fontsize=12)
    ax.set_ylabel('Scattering rate (1/ns)', fontsize=12)
    legend = ["Sidewalls diffuse", "Sidewalls specular", "Top & bottom diffuse", "Top & bottom specular",
              "Holes diffuse", "Holes specular", "Hot side", "Internal", "Pillars diffuse", "Pillars specular"]
    ax.legend(legend, loc='upper right')
    plt.yscale('log')
    ax.set_ylim(bottom=1.0)
    fig.savefig("Scattering rates.pdf", dpi=300, format='pdf', bbox_inches="tight")
    if cf.plots_in_terminal: plt.show()

    # Save the file:
    filename = "Data/Scattering rates.csv"
    header1 = "Y [um], Sidewalls diffuse [1/ns], Sidewalls specular [1/ns], Top & bottom diffuse [1/ns], "
    header2 = "Top & bottom specular [1/ns], Holes diffuse [1/ns], Holes specular [1/ns], Hot side [1/ns], "
    header3 = "Internal [1/ns], Pillars diffuse [1/ns], Pillars specular [1/ns]"
    header = header1 + header2 + header3
    data = np.vstack((segments, all_scattering_rates)).T
    np.savetxt(filename, data, fmt='%1.2e', delimiter=",", header=header)


def plot_data():
    """Create plots of various distributions"""
    plot_angle_distribution()
    plot_free_path_distribution()
    plot_frequency_distribution()
    plot_wavelength_distribution()
    plot_travel_time_distribution()
    plot_detected_frequency_distribution()
    plot_velocity_distribution()
    plot_time_in_segments()
    plot_thermal_conductivity()
    plot_temperature_profile()
    plot_heat_flux_profile()
    plot_thermal_map()
    plot_trajectories()
    plot_scattering_statistics()
    if cf.output_scattering_map:
        plot_scattering_map()
