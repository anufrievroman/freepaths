"""Main module of the program"""

import os
import sys
import time
import shutil

# Modules:
from parameters import *
from output import *
from analysis import plot_data
from run_phonon import run_phonon
from phonon import Phonon
from flight import Flight
from data import ScatteringData, GeneralData, SegmentData, PathData
from progress import Progress
from materials import Material
from maps import ScatteringMap, ThermalMaps


def main():
    """This is the main function, which works under Debye approximation.
    It should be used to simulate phonon paths at low temperatures"""

    print(f'Simulation for {OUTPUT_FOLDER_NAME} started.')
    start_time = time.time()
    progress = Progress()

    # Initiate data structures:
    material = Material(MEDIA)
    scatter_stats = ScatteringData()
    general_stats = GeneralData()
    segment_stats = SegmentData()
    path_stats = PathData()
    scatter_maps = ScatteringMap()
    thermal_maps = ThermalMaps()

    # For each phonon:
    for index in range(NUMBER_OF_PHONONS):
        progress.render(index)

        # Initiate a phonon and its flight:
        phonon = Phonon(material)
        flight = Flight(phonon)

        # Run this phonon through the structure:
        run_phonon(phonon, flight, scatter_stats, segment_stats, thermal_maps, scatter_maps, material)

        # Record the properties returned for this phonon:
        general_stats.save_phonon_data(phonon)
        general_stats.save_flight_data(flight)

        # Record trajectories of the first N phonons:
        if index < OUTPUT_TRAJECTORIES_OF_FIRST:
            path_stats.save_phonon_path(flight)

    # Run additional calculations:
    thermal_maps.calculate_thermal_conductivity()

    # Create the folder if it does not exist and copy parameters.py there:
    if not os.path.exists("Results/" + OUTPUT_FOLDER_NAME):
        os.makedirs("Results/" + OUTPUT_FOLDER_NAME)
        os.makedirs("Results/" + OUTPUT_FOLDER_NAME + '/Data')
    shutil.copy('parameters.py', "Results/" + OUTPUT_FOLDER_NAME)
    os.chdir("Results/" + OUTPUT_FOLDER_NAME)

    # Save data into files:
    general_stats.write_into_files()
    scatter_stats.write_into_files()
    segment_stats.write_into_files()
    path_stats.write_into_files()
    thermal_maps.write_into_files()
    scatter_maps.write_into_files()

    # Output general information:
    output_general_information(start_time)
    output_scattering_information(scatter_stats)

    # Analyze and plot the data:
    sys.stdout.write("\rAnalyzing the data...")
    plot_data()

    sys.stdout.write(f'\rSee the results in "Results/{OUTPUT_FOLDER_NAME}" folder.\n')
    sys.stdout.write("\rThank you for using FreePATHS.\n")


if __name__ == "__main__":
    main()
