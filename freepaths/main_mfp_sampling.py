"""Module to run thermal conductivity calculations by integrating the phonon dispersion"""

import os
import sys
import time
import shutil
import scipy
import math

# Modules:
from freepaths.animation import create_animation
from freepaths.config import cf
from freepaths.run_phonon import run_phonon
from freepaths.phonon import Phonon, Polarization
from freepaths.flight import Flight
from freepaths.data import ScatteringData, GeneralData, SegmentData, PathData
from freepaths.progress import Progress
from freepaths.materials import Material
from freepaths.maps import ScatteringMap, ThermalMaps
from freepaths.output_info import output_general_information, output_scattering_information
from freepaths.output_plots import plot_data


def main(input_file):
    """This is the main function, which integrates phonon dispersion to get thermal conductivity"""

    print(f'Mean free path sampling of {cf.output_folder_name}')
    start_time = time.time()
    progress = Progress()

    # Initiate data structures:
    material = Material(cf.media, num_points=cf.number_of_phonons+1)
    scatter_stats = ScatteringData()
    general_stats = GeneralData()
    segment_stats = SegmentData()
    path_stats = PathData()
    scatter_maps = ScatteringMap()
    thermal_maps = ThermalMaps()

    thermal_conductivity = 0

    # For each polarization branch:
    for polarization in [Polarization.LA, Polarization.TA, Polarization.TA]:
        sys.stdout.write(f"\rIntegrating {polarization.name} branch.\n")

        # For each phonon:
        for index in range(cf.number_of_phonons):

            # Wave vector:
            k_vector = (material.dispersion[index+1, 0] + material.dispersion[index, 0]) / 2
            d_k_vector = (material.dispersion[index+1, 0] - material.dispersion[index, 0])

            # Initiate a phonon and its flight:
            phonon = Phonon(material, polarization, index)
            flight = Flight(phonon)

            # Run this phonon through the structure:
            run_phonon(phonon, flight, scatter_stats, segment_stats, thermal_maps, scatter_maps, material)

            # Record the properties returned for this phonon:
            general_stats.save_phonon_data(phonon)
            general_stats.save_flight_data(flight)

            # Record trajectories of the first N phonons:
            if index < cf.output_trajectories_of_first:
                path_stats.save_phonon_path(flight)

            # Heat capacity, Ref. PRB 88 155318 (2013):
            omega = 2 * math.pi * phonon.f
            part = scipy.constants.hbar * omega / (scipy.constants.k * cf.temp)
            c_p = scipy.constants.k * part**2 * math.exp(part) / (math.exp(part) - 1)**2

            # Thermal conductivity, Ref. Phys. Rev. 132 2461 (1963):
            relax_time = flight.mean_free_path/phonon.speed
            thermal_conductivity += (1/(6*(math.pi**2)))*c_p*(phonon.speed**2)*relax_time*(k_vector**2)*d_k_vector

    # Run additional calculations:
    thermal_maps.calculate_thermal_conductivity()

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

    # Output general information:
    output_general_information(start_time)
    output_scattering_information(scatter_stats)

    # Analyze and plot the data:
    sys.stdout.write("\rAnalyzing the data...")
    plot_data()

    # Generate animation of phonon paths:
    if cf.output_path_animation:
        sys.stdout.write("\rGenerating path animation...")
        create_animation()

    sys.stdout.write(f'\rSee the results in "Results/{cf.output_folder_name}" folder.\n')
    sys.stdout.write(f"\rThermal conductivity = {thermal_conductivity}\n")
    sys.stdout.write("\rThank you for using FreePATHS.\n")
