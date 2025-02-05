"""Module that calculates and outputs vaious plots and distributions from the saved files"""

import logging
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LogNorm
from matplotlib import font_manager

from freepaths.config import cf
from freepaths.output_structure import draw_structure_top_view, draw_structure_side_view
from freepaths.materials import Si, SiC, Graphite
import matplotlib.pyplot as plt

# Style of the plots:
all_fonts = font_manager.get_font_names()
if 'Arial' in all_fonts:
    plt.rcParams['font.family'] = 'Arial'
else:
    # Use a generic font or specify multiple options
    logging.warning("Arial font not available. Falling back on default sans-serif font")
    plt.rcParams['font.family'] = ['sans-serif']
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.titlesize'] = 8
plt.rcParams['lines.linewidth'] = 1.0
plt.rcParams['xtick.direction'] = "in"
plt.rcParams['ytick.direction'] = "in"
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['xtick.minor.visible'] = False
plt.rcParams['ytick.minor.visible'] = False
plt.rcParams['legend.frameon'] = False
plt.rcParams['figure.autolayout'] = True
plt.rcParams['figure.figsize'] = [5, 3.5]
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['grid.linewidth'] = 0.5


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
    hole_diff_angles = np.loadtxt("Data/All hole diffuse scattering angles.csv", dtype='float', encoding='utf-8')
    hole_spec_angles = np.loadtxt("Data/All hole specular scattering angles.csv", dtype='float', encoding='utf-8')
    distribution = np.zeros((360, 3))
    distribution[:, 0] = range(-180, 180)
    exit_angles = all_exit_angles[all_exit_angles != 0]
    distribution[:, 1], _ = np.histogram(np.degrees(exit_angles), 360, range=(-180, 180))
    distribution[:, 2], _ = np.histogram(np.degrees(initial_angles), 360, range=(-180, 180))
    return distribution


def scattering_angle_distribution_calculation():
    """Analyse scattering phonon angles and create their distribution"""
    hole_diff_angles = np.loadtxt("Data/All hole diffuse scattering angles.csv", dtype='float', encoding='utf-8')
    hole_spec_angles = np.loadtxt("Data/All hole specular scattering angles.csv", dtype='float', encoding='utf-8')
    distribution = np.zeros((360, 3))
    distribution[:, 0] = range(-180, 180)
    distribution[:, 1], _ = np.histogram(np.degrees(hole_diff_angles), 360, range=(-180, 180))
    distribution[:, 2], _ = np.histogram(np.degrees(hole_spec_angles), 360, range=(-180, 180))
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


def cumulative_conductivity_calculation():
    """Calculate cumulative function of thermal conductivity vs mean freee path"""
    mfp = np.loadtxt("Data/All mean free paths.csv", encoding='utf-8')
    kappa = np.loadtxt("Data/All thermal conductivities.csv", encoding='utf-8')
    sorted_indices = np.argsort(mfp)
    sorted_mfp = mfp[sorted_indices]
    sorted_kappa = kappa[sorted_indices]
    cumulative_thermal_conductivity = np.cumsum(sorted_kappa)
    return sorted_mfp, cumulative_thermal_conductivity


def plot_cumulative_thermal_conductivity(mfp_sampling):
    """Plot distribution cumulative thermal conductivity vs mean free path"""
    if not mfp_sampling:
        return
    mfp, kappa = cumulative_conductivity_calculation()
    fig, ax = plt.subplots()
    ax.plot(mfp * 1e6, kappa, 'royalblue')
    ax.set_xlabel('Mean free path (μm)')
    ax.set_ylabel('Cumulative thermal conductivity (W/m·K)')
    ax.set_xscale('log')
    fig.savefig("Distribution of thermal conductivity.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)
    np.savetxt('Data/Distribution of thermal conductivity.csv', np.vstack((mfp, kappa)).T, fmt='%1.3e', delimiter=",")


def plot_angle_distribution():
    """Plot distribution of initial and exit angles"""
    angle_distributions = angle_distribution_calculation()
    fig, ax = plt.subplots()
    ax.plot(angle_distributions[:, 0], angle_distributions[:, 1], 'royalblue')
    ax.plot(angle_distributions[:, 0], angle_distributions[:, 2], 'deeppink')
    ax.set_xlabel('Angle (degree)')
    ax.set_ylabel('Number of phonons')
    ax.set_ylim(bottom=0)
    ax.legend(["At cold side", "At hot side"])
    fig.savefig("Distribution of angles.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)
    np.savetxt('Data/Distribution of angles.csv', angle_distributions, fmt='%1.3e', delimiter=",")


def plot_scattering_angle_distribution():
    """Plot distribution of hole scattering angles"""
    angle_distributions = scattering_angle_distribution_calculation()
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.grid(zorder=0)
    ax.plot(np.deg2rad(angle_distributions[:, 0]), angle_distributions[:, 1], 'royalblue', label="Diffuse", zorder=2)
    ax.fill_between(np.deg2rad(angle_distributions[:, 0]), angle_distributions[:, 1], 0, alpha=0.2, zorder=2)
    ax.plot(np.deg2rad(angle_distributions[:, 0]), angle_distributions[:, 2], 'deeppink', label="Specular", zorder=3)
    ax.fill_between(np.deg2rad(angle_distributions[:, 0]), angle_distributions[:, 2], 0, alpha=0.2, zorder=3)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.legend(facecolor='white', framealpha=1, ncols=2, loc="lower center", bbox_to_anchor=(0.5, -0.17))
    ax.set_title('Number of scattered phonons per angle')
    fig.savefig("Distribution of hole scattering angles.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)
    np.savetxt('Data/Distribution of hole scattering angles.csv', angle_distributions, fmt='%1.3e', delimiter=",")


def plot_free_path_distribution():
    """Plot distribution of free path"""
    filename = "Data/All free paths.csv"
    free_path_distribution = distribution_calculation(filename, None, cf.number_of_nodes)
    fig, ax = plt.subplots()
    ax.plot(free_path_distribution[:, 0] * 1e6, free_path_distribution[:, 1], 'royalblue')
    ax.set_xscale('log')
    ax.set_xlabel('Free flights (μm)')
    ax.set_ylabel('Number of flights')
    fig.savefig("Distribution of free paths.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)
    np.savetxt('Data/Distribution of free paths.csv', free_path_distribution, fmt='%1.3e', delimiter=",")


def plot_frequency_distribution():
    """Plot distribution of frequencies"""
    filename = "Data/All initial frequencies.csv"
    frequency_distribution = distribution_calculation(filename, None, cf.number_of_nodes)
    fig, ax = plt.subplots()
    ax.plot(frequency_distribution[:, 0] * 1e-12, frequency_distribution[:, 1], 'royalblue')
    ax.set_xlabel('Frequency (THz)')
    ax.set_ylabel('Number of phonons')
    fig.savefig("Distribution of initial frequencies.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)
    np.savetxt('Data/Distribution of initial frequencies.csv', frequency_distribution, fmt='%1.3e', delimiter=",")


def plot_wavelength_distribution():
    """Plot distribution of wavelength"""
    wavelength_distribution = wavelength_distribution_calculation(cf.number_of_nodes)
    fig, ax = plt.subplots()
    ax.plot(wavelength_distribution[:, 0] * 1e9, wavelength_distribution[:, 1], 'royalblue')
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Number of phonons')
    fig.savefig("Distribution of wavelengths.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)
    np.savetxt('Data/Distribution of wavelengths.csv', wavelength_distribution, fmt='%1.3e', delimiter=",")


def plot_travel_time_distribution():
    """Plot distribution of wavelength"""
    travel_time_distribution = distribution_calculation("Data/All travel times.csv", None, cf.number_of_nodes)
    fig, ax = plt.subplots()
    ax.plot(travel_time_distribution[:, 0] * 1e9, travel_time_distribution[:, 1], 'royalblue')
    ax.set_xscale('log')
    ax.set_xlabel('Travel time (ns)')
    ax.set_ylabel('Number of phonons')
    fig.savefig("Distribution of travel times.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)
    np.savetxt('Data/Distribution of travel times.csv', travel_time_distribution, fmt='%1.3e', delimiter=",")


def plot_mean_free_path_distribution():
    """Plot distribution of MFP per phonon"""
    mean_free_path_distribution = distribution_calculation("Data/All mean free paths.csv", None, cf.number_of_nodes)
    fig, ax = plt.subplots()
    ax.plot(mean_free_path_distribution[:, 0] * 1e9, mean_free_path_distribution[:, 1], 'royalblue')
    ax.set_xscale('log')
    ax.set_xlabel('Mean free path (nm)')
    ax.set_ylabel('Number of phonons')
    fig.savefig("Distribution of MFPs.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)
    np.savetxt('Data/Distribution of MFPs.csv', mean_free_path_distribution, fmt='%1.3e', delimiter=",")


def plot_velocity_distribution():
    """Plot distribution of group velocities"""
    fig, ax = plt.subplots()
    speeds = np.loadtxt("Data/All group velocities.csv")
    frequencies = np.loadtxt("Data/All initial frequencies.csv")
    ax.plot(frequencies, speeds, '.', markersize=2, c='royalblue')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Group velocity (m/s)')
    fig.savefig('Group velocities.pdf', format='pdf', bbox_inches="tight")
    plt.close(fig)


def plot_time_in_segments():
    """Plot time spent in segments"""
    fig, ax = plt.subplots()
    segment, time = np.genfromtxt("Data/Time spent in segments.csv", unpack=True, delimiter=',', usecols=(0, 1), skip_header=1)
    ax.plot(segment, time, '-o', markersize=2, c='royalblue')
    ax.set_xlabel('Y (μm)')
    ax.set_ylabel('Time spent (ns)')
    fig.savefig("Time spent in segments.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)


def plot_thermal_conductivity():
    """Plot thermal conductivity against time segment"""
    fig, ax = plt.subplots()
    time, eff_tc, mat_tc = np.genfromtxt("Data/Thermal conductivity.csv", unpack=True, delimiter=',', usecols=(0, 1, 2), skip_header=1)
    av_eff_tc, av_mat_tc, std_eff_tc, std_mat_tc, av_start, av_end = np.genfromtxt("Data/Average thermal conductivity.csv", unpack=True, delimiter=',', skip_header=1)

    # Plot thermal conductivity vs time:
    ax.plot(time, mat_tc, '-o', markersize=2, linewidth=1, color='deeppink', label='κ$_{mat}$')
    ax.plot(time, eff_tc, '-o', markersize=2, linewidth=1, color='royalblue', label='κ$_{eff}$')

    # Add labels:
    ax.text(av_start/2, max(mat_tc), 'Stabilization period', color='grey', fontsize=8, ha='center')
    ax.text(av_start + (av_end-av_start)/2, max(mat_tc), 'Measurement period', color='grey', fontsize=8, ha='center')

    # Plot the range of averaging:
    ax.axvline(x=av_start, color='grey', linestyle='--')
    ax.axvline(x=av_end, color='grey', linestyle='--')

    # Plot the standard deviation range:
    rectangle_mat = patches.Rectangle((av_start, av_mat_tc-std_mat_tc), av_end-av_start, 2*std_mat_tc, linewidth=0, facecolor='deeppink', alpha=0.2)
    rectangle_eff = patches.Rectangle((av_start, av_eff_tc-std_eff_tc), av_end-av_start, 2*std_eff_tc, linewidth=0, facecolor='royalblue', alpha=0.2)
    ax.add_patch(rectangle_eff)
    ax.add_patch(rectangle_mat)

    # Plot the averaged values in the averaging range:
    ax.plot([av_start, av_end], [av_mat_tc, av_mat_tc], '-', markersize=2, linewidth=2, color='deeppink', label=r'$\overline{\kappa}_{mat}$')
    ax.plot([av_start, av_end], [av_eff_tc, av_eff_tc], '-', markersize=2, linewidth=2, color='royalblue', label=r'$\overline{\kappa}_{eff}$')

    ax.set_ylabel('Thermal conductivity (W/m·K)')
    ax.set_xlabel('Time (ns)')

    ax.set_xlim(left=0.0)
    ax.legend()
    kappa_mat_str = r"$\overline{\kappa}_{mat}$"
    kappa_eff_str = r"$\overline{\kappa}_{eff}$"
    ax.set_title(f'{kappa_mat_str} = {av_mat_tc} ± {std_mat_tc} W/m·K, {kappa_eff_str} = {av_eff_tc} ± {std_eff_tc} W/m·K', color='grey')
    fig.savefig("Thermal conductivity.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)


def plot_temperature_profile():
    """Plot profile of temperature for each time segment"""
    fig, ax = plt.subplots()
    data = np.genfromtxt("Data/Temperature profiles y.csv", unpack=True, delimiter=',', skip_header=1, encoding='utf-8')
    for timeframe_num in range(len(data) - 1):
        ax.plot(data[0][1:], data[timeframe_num + 1][1:], linewidth=1, label=f'Time frame {timeframe_num+1}')
    ax.set_xlabel('Y (μm)')
    ax.set_ylabel('Temperature (K)')
    ax.legend()
    fig.savefig("Temperature profile.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)


def plot_heat_flux_profile():
    """Plot profile of heat flux for each time segment"""

    # Effective heat flux:
    fig, ax = plt.subplots()
    data = np.genfromtxt("Data/Heat flux profiles y.csv", unpack=True, delimiter=',', skip_header=1, encoding='utf-8')
    for timeframe_num in range((len(data) - 1) // 2):
        ax.plot(data[0], data[timeframe_num + 1], linewidth=1, label=f'Time frame {timeframe_num+1}')
    ax.set_xlabel('Y (μm)')
    ax.set_ylabel('Heat flux (W/m²)')
    ax.legend()
    fig.savefig("Heat flux profile effective.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)

    # Material heat flux:
    fig, ax = plt.subplots()
    for timeframe_num in range((len(data) - 1) // 2):
        ax.plot(data[0], data[2*timeframe_num + 1], linewidth=1, label=f'Time frame {timeframe_num+1}')
    ax.set_xlabel('Y (μm)')
    ax.set_ylabel('Heat flux (W/m²)')
    ax.legend()
    fig.savefig("Heat flux profile material.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)

def plot_thermal_map():
    """Plot thermal map as color map"""
    fig = plt.figure()
    thermal_map = np.genfromtxt("Data/Thermal map.csv", unpack=False, delimiter=',', skip_header=0, encoding='utf-8')
    thermal_map = np.flipud(thermal_map)
    minimum_of_colorbar = 1e-20  # Cannot be zero!
    boundaries = [(-cf.width / 2) * 1e6, (cf.width / 2) * 1e6, 0, cf.length * 1e6]
    plt.imshow(thermal_map, cmap='hot', interpolation='none', extent=boundaries,
               norm=LogNorm(vmin=minimum_of_colorbar, vmax=np.amax(thermal_map)))
    plt.xlabel('x (μm)')
    plt.ylabel('y (μm)')
    cbar = plt.colorbar(shrink=0.6)
    cbar.set_label('Energy density', rotation=90)
    fig.savefig("Thermal map.pdf", bbox_inches="tight")
    plt.close(fig)

def plot_pixel_volumes():
    """Plot the pixel volumes as 2D map"""
    fig = plt.figure()
    pixel_volumes = np.genfromtxt("Data/Pixel volumes.csv", unpack=False, delimiter=',', skip_header=0, encoding='utf-8')
    pixel_volumes = np.flipud(pixel_volumes)
    boundaries = [(-cf.width / 2) * 1e6, (cf.width / 2) * 1e6, 0, cf.length * 1e6]
    plt.imshow(pixel_volumes, cmap='hot', interpolation='none', extent=boundaries)
    plt.xlabel('x (μm)')
    plt.ylabel('y (μm)')
    fig.savefig("Pixel volumes.pdf", bbox_inches="tight")
    plt.close(fig)

def plot_heat_flux_map(file, label, units="a.u."):
    """Plot heat flux map as color map"""
    fig = plt.figure()
    heat_flux_map = np.genfromtxt(file, unpack=False, delimiter=',', skip_header=0, encoding='utf-8')
    heat_flux_map = np.flipud(heat_flux_map)
    minimum_of_colorbar = 1e5
    boundaries = [(-cf.width / 2) * 1e6, (cf.width / 2) * 1e6, 0, cf.length * 1e6]
    plt.imshow(heat_flux_map, cmap='jet', interpolation='none', extent=boundaries,
               norm=LogNorm(vmin=minimum_of_colorbar, vmax=np.amax(heat_flux_map)))
    plt.xlabel('x (μm)', fontsize=10)
    plt.ylabel('y (μm)', fontsize=10)
    cbar = plt.colorbar(shrink=0.6)
    cbar.set_label(f"{label} ({units})", rotation=90, fontsize=10)
    fig.savefig(f"{label}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_scattering_map():
    """Plot the map of scattering events"""
    if not cf.output_scattering_map:
        return
    fig, ax = plt.subplots()
    filename = "Data/Scattering map.csv"
    spec_x, spec_y, diff_x, diff_y, int_x, int_y = np.genfromtxt(filename, unpack=True, delimiter=',', skip_header=1, encoding='utf-8')
    ax.plot(np.trim_zeros(diff_x, 'b')*1e6, np.trim_zeros(diff_y, 'b')*1e6, 'o', color='b', markersize=0.1, alpha=0.3)
    ax.plot(np.trim_zeros(spec_x, 'b')*1e6, np.trim_zeros(spec_y, 'b')*1e6, 'o', color='g', markersize=0.1, alpha=0.3)
    ax.plot(np.trim_zeros(int_x, 'b')*1e6, np.trim_zeros(int_y, 'b')*1e6, 'o', color='r', markersize=0.1, alpha=0.3)
    ax.set_xlabel('x (μm)')
    ax.set_ylabel('y (μm)')
    ax.legend(["Diffuse", "Specular", "Internal"])
    ax.set_aspect('equal', 'datalim')
    fig.savefig("Scattering map.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_trajectories():
    """Plot the phonon trajectories"""

    data = np.genfromtxt("Data/Phonon paths.csv", unpack=False, delimiter=',', skip_header=1, encoding='utf-8')

    # Create XY plot:
    fig, ax = plt.subplots()

    # Draw structure:
    patches = draw_structure_top_view(cf, color_holes='white', color_back=cf.output_structure_color)
    for patch in patches:
        ax.add_patch(patch)

    # Draw paths:
    for index in range(cf.output_trajectories_of_first):
        x_coordinates = np.trim_zeros(data[:, 3 * index], trim='b')
        y_coordinates = np.trim_zeros(data[:, 3 * index + 1], trim='b')
        max_steps = min([x_coordinates.shape[0], y_coordinates.shape[0]])
        ax.plot(x_coordinates[:max_steps], y_coordinates[:max_steps], linewidth=0.2)

    # Set labels:
    ax.set_xlabel('X (μm)')
    ax.set_ylabel('Y (μm)')
    ax.set_aspect('equal', 'datalim')
    fig.savefig("Phonon paths XY.pdf", dpi=600, format='pdf', bbox_inches="tight")
    plt.close(fig)

    # Create YZ plot:
    fig, ax = plt.subplots()

    # Draw structure:
    patches = draw_structure_side_view(cf, color_holes='white', color_back=cf.output_structure_color)
    for patch in patches:
        ax.add_patch(patch)

    for index in range(cf.output_trajectories_of_first):
        y_coordinates = np.trim_zeros(data[:, 3 * index + 1], trim='b')
        num_of_points = len(y_coordinates)
        z_coordinates = data[:num_of_points, 3 * index + 2]
        max_steps = min([y_coordinates.shape[0], z_coordinates.shape[0]])
        ax.plot(y_coordinates[:max_steps], z_coordinates[:max_steps], linewidth=0.2)
    ax.set_xlabel('Y (μm)')
    ax.set_ylabel('Z (μm)')
    ax.set_aspect('equal', 'datalim')
    fig.savefig("Phonon paths YZ.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)


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
        if scattering_type == 6: # Skip hot side scattering
            continue
        # Calculate and plot scattering events per second:
        scattering_rate = [events / time for events, time in zip(scattering_data[scattering_type, :], time_spent)]
        all_scattering_rates.append(scattering_rate)
        ax.plot(segments, scattering_rate,  '-o', markersize=2)

    ax.set_xlabel('Y (μm)')
    ax.set_ylabel('Scattering rate (1/ns)')
    legend = ["Sidewalls diffuse", "Sidewalls specular", "Top & bottom diffuse", "Top & bottom specular",
              "Holes diffuse", "Holes specular", "Internal", "Pillars diffuse", "Pillars specular"]
    # ax.legend(legend, loc='upper right')
    ax.legend(legend, loc='lower center', ncol=2, fancybox=True)
    plt.yscale('log')
    ax.set_ylim(bottom=1.0)
    fig.savefig("Scattering rates.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)

    # Save the file:
    filename = "Data/Scattering rates.csv"
    header1 = "Y [um], Sidewalls diffuse [1/ns], Sidewalls specular [1/ns], Top & bottom diffuse [1/ns], "
    header2 = "Top & bottom specular [1/ns], Holes diffuse [1/ns], Holes specular [1/ns], Hot side [1/ns], "
    header3 = "Internal [1/ns], Pillars diffuse [1/ns], Pillars specular [1/ns]"
    header = header1 + header2 + header3
    data = np.vstack((segments, all_scattering_rates)).T
    np.savetxt(filename, data, fmt='%1.2e', delimiter=",", header=header)


def plot_structure():
    """Plot the structure with all the elements"""

    # Create XY plot:
    fig, ax = plt.subplots()

    # Draw structures:
    patches = draw_structure_top_view(cf, color_holes='black', color_back='royalblue')
    for patch in patches:
        ax.add_patch(patch)

    # Set labels:
    ax.set_xlabel('X (μm)')
    ax.set_ylabel('Y (μm)')
    ax.set_aspect('equal', 'datalim')
    ax.set_xlim(1.2*-cf.width*1e6, 1.2*cf.width*1e6)
    ax.set_ylim(0, cf.length*1e6)
    fig.savefig("Structure XY.pdf", dpi=600, format='pdf', bbox_inches="tight")
    plt.close(fig)


def plot_material_properties():
    """Plot phonon dispersion and display some other material properties"""

    # Initialize the material:
    if cf.media == "Si":
        material = Si(cf.temp)
    elif cf.media == "SiC":
        material = SiC(cf.temp)
    elif cf.media == "Graphite":
        material = Graphite(cf.temp)
    else:
        logging.error(f"Material {cf.media} is not supported")
        sys.exit()

    # Plot phonon dispersion:
    fig, ax = plt.subplots()
    for index in range(material.dispersion.shape[1] - 1):
        ax.plot(material.dispersion[:,0], material.dispersion[:,index + 1] * 1e-12, linewidth=1, label=f'{index}')
    ax.set_xlabel('Wavevector (1/m)')
    ax.set_ylabel('Frequency (THz)')
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0)

    # Add material properties:
    ax.set_title(f'{material.name},  T = {cf.temp} K,  C$_p$ = {material.heat_capacity:.3f} J/kg·K,  ρ = {material.density} kg/m³', color="grey")

    fig.savefig("Material properties.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)


def plot_data(mfp_sampling=False):
    """Create plots of various distributions, maps, profiles, and other quantities"""
    function_list = [
        plot_structure,
        plot_trajectories,
        plot_angle_distribution,
        plot_scattering_angle_distribution,
        plot_free_path_distribution,
        plot_frequency_distribution,
        plot_wavelength_distribution,
        plot_travel_time_distribution,
        plot_mean_free_path_distribution,
        plot_velocity_distribution,
        plot_time_in_segments,
        plot_thermal_conductivity,
        plot_temperature_profile,
        plot_heat_flux_profile,
        plot_thermal_map,
        plot_pixel_volumes,
        plot_scattering_statistics,
        plot_scattering_map,
        plot_material_properties,
    ]

    # Run main functions and handle exceptions:
    for func in function_list:
        try:
            func()
        except Exception as e:
            logging.warning(f"Function {func.__name__} failed: {e}")

    # Run additional functions:
    plot_cumulative_thermal_conductivity(mfp_sampling)
    plot_heat_flux_map(file="Data/Heat flux map xy.csv", label="Heat flux map", units="W/m²")
    plot_heat_flux_map(file="Data/Heat flux map x.csv", label="Heat flux map x", units="W/m²")
    plot_heat_flux_map(file="Data/Heat flux map y.csv", label="Heat flux map y", units="W/m²")
