"""Module to run thermal conductivity calculations by integrating the phonon dispersion"""

import os
import sys
import time
import shutil
import scipy
import math
import logging
from colorama import Fore, Style

# Modules:
from freepaths.animation import create_animation
from freepaths.config import cf
from freepaths.run_phonon import run_phonon
from freepaths.phonon import Phonon
from freepaths.flight import Flight
from freepaths.data import ScatteringData, GeneralData, SegmentData, PathData, TriangleScatteringData
from freepaths.progress import Progress
from freepaths.materials import Si, SiC, Graphite
from freepaths.maps import ScatteringMap, ThermalMaps
from freepaths.output_info import output_general_information, output_scattering_information, output_parameter_warnings
from freepaths.output_plots import plot_data


def main(input_file):
    """This is the main function, which integrates phonon dispersion to get thermal conductivity"""

    print(f'Mean free path sampling of {Fore.GREEN}{cf.output_folder_name}{Style.RESET_ALL}')
    start_time = time.time()
    progress = Progress()

    # Initialize the material:
    if cf.media == "Si":
        material = Si(cf.temp, num_points=cf.number_of_phonons +1)
    elif cf.media == "SiC":
        material = SiC(cf.temp, num_points=cf.number_of_phonons+1)
    elif cf.media == "Graphite":
        material = Graphite(cf.temp, num_points=cf.number_of_phonons+1)
    else:
        logging.error(f"Material {cf.media} is not supported")
        sys.exit()

    # Initiate data structures:
    scatter_stats = ScatteringData()
    general_stats = GeneralData()
    segment_stats = SegmentData()
    places_stats = TriangleScatteringData()
    path_stats = PathData()
    scatter_maps = ScatteringMap()
    thermal_maps = ThermalMaps()

    total_thermal_conductivity = 0.0

    # For each polarization branch:
    for branch_number in range(3):
        sys.stdout.write(f"\rIntegrating branch number {branch_number+1}.\n")

        # For each phonon:
        for index in range(cf.number_of_phonons):

            # Wave vector:
            k_vector = (material.dispersion[index+1, 0] + material.dispersion[index, 0]) / 2
            d_k_vector = (material.dispersion[index+1, 0] - material.dispersion[index, 0])

            # Initiate a phonon and its flight:
            phonon = Phonon(material, branch_number, index)
            flight = Flight(phonon)

            # Run this phonon through the structure:
            run_phonon(phonon, flight, scatter_stats, places_stats, segment_stats, thermal_maps, scatter_maps, material)

            # Heat capacity, Ref. PRB 88 155318 (2013):
            omega = 2 * math.pi * phonon.f
            part = scipy.constants.hbar * omega / (scipy.constants.k * cf.temp)
            c_p = scipy.constants.k * part**2 * math.exp(part) / (math.exp(part) - 1)**2

            # Thermal conductivity, Ref. Phys. Rev. 132 2461 (1963):
            mean_relax_time = flight.mean_free_path/phonon.speed
            flight.thermal_conductivity = (1/(6*(math.pi**2)))*c_p*(phonon.speed**2)*mean_relax_time*(k_vector**2)*d_k_vector
            total_thermal_conductivity += flight.thermal_conductivity

            # Record the properties returned for this phonon:
            general_stats.save_phonon_data(phonon)
            general_stats.save_flight_data(flight)

            # Record trajectories of the first N phonons:
            if index < cf.output_trajectories_of_first:
                path_stats.save_phonon_path(flight)

    # Run additional calculations:
    thermal_maps.calculate_thermal_conductivity()
    thermal_maps.calculate_weighted_flux()
    thermal_maps.calculate_heat_flux_modulus()

    # Create the folder if it does not exist and copy input file there:
    if not os.path.exists("Results/" + cf.output_folder_name):
        os.makedirs("Results/" + cf.output_folder_name)
        os.makedirs("Results/" + cf.output_folder_name + '/Data')
    if input_file:
        shutil.copy(input_file, "Results/" + cf.output_folder_name)
    os.chdir("Results/" + cf.output_folder_name)

    # Save data into files:
    general_stats.write_into_files()
    scatter_stats.write_into_files()
    segment_stats.write_into_files()
    thermal_maps.write_into_files()
    scatter_maps.write_into_files()
    path_stats.write_into_files()

    # Generate animation of phonon paths:
    if cf.output_path_animation:
        create_animation()

    # Analyze and plot the data:
    sys.stdout.write("\rAnalyzing the data...")
    plot_data(mfp_sampling=True)

    # Output general information:
    output_general_information(start_time)
    output_scattering_information(scatter_stats)
    output_parameter_warnings()

    sys.stdout.write(f'\rSee the results in {Fore.GREEN}Results/{cf.output_folder_name}{Style.RESET_ALL}\n')
    sys.stdout.write(f"\rThermal conductivity = {Fore.GREEN}{total_thermal_conductivity:.5f}{Style.RESET_ALL} W/m·K\n")
    sys.stdout.write(f"\r{Fore.BLUE}Thank you for using FreePATHS{Style.RESET_ALL}\n\n")
