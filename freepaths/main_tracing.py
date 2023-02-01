"""Main module of the program"""

import os
import sys
import time
import shutil

# Modules:
from freepaths.config import cf
from freepaths.run_phonon import run_phonon
from freepaths.phonon import Phonon
from freepaths.flight import Flight
from freepaths.data import ScatteringData, GeneralData, SegmentData, PathData
from freepaths.progress import Progress
from freepaths.materials import Material
from freepaths.maps import ScatteringMap, ThermalMaps
from freepaths.output_info import output_general_information, output_scattering_information
from freepaths.animation import create_animation
from freepaths.output_plots import plot_data


def main(input_file):
    """This is the main function, which works under Debye approximation.
    It should be used to simulate phonon paths at low temperatures"""

    print(f'Simulation of {cf.output_folder_name}')
    start_time = time.time()
    progress = Progress()

    # Initiate data structures:
    material = Material(cf.media)
    scatter_stats = ScatteringData()
    general_stats = GeneralData()
    segment_stats = SegmentData()
    path_stats = PathData()
    scatter_maps = ScatteringMap()
    thermal_maps = ThermalMaps()

    # For each phonon
    for index in range(cf.number_of_phonons):
        progress.render(index, cf.number_of_phonons)

        # Initiate a phonon and its flight:
        phonon = Phonon(material)
        flight = Flight(phonon)

        # Run this phonon through the structure:
        run_phonon(phonon, flight, scatter_stats, segment_stats, thermal_maps, scatter_maps, material)

        # Record the properties returned for this phonon:
        general_stats.save_phonon_data(phonon)
        general_stats.save_flight_data(flight)

        # Record trajectories of the first N phonons:
        if index < cf.output_trajectories_of_first:
            path_stats.save_phonon_path(flight)

    # Run additional calculations:
    thermal_maps.calculate_thermal_conductivity()

    # Create the folder if it does not exist and copy input file there:
    if not os.path.exists(f"Results/{cf.output_folder_name}"):
        os.makedirs(f"Results/{cf.output_folder_name}")
        os.makedirs(f"Results/{cf.output_folder_name}/Data")
    if cf.output_path_animation and not os.path.exists(f"Results/{cf.output_folder_name}/Frames"):
        os.makedirs(f"Results/{cf.output_folder_name}/Frames")
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
    sys.stdout.write("\rThank you for using FreePATHS.\n")
