"""Main module of the program"""

import os
import sys
import time
import shutil

# Modules:
from parameters import *
from output import *
from scattering import *
from analysis import plot_data
from phonon import Phonon
from flight import Flight
from data import ScatteringData, GeneralData, SegmentData, PathData
from progress import Progress
from materials import Material
from maps import ScatteringMap, ThermalMaps
from events import ScatteringTypes


def run_phonon(phonon, flight, scatter_stats, segment_stats, thermal_maps, scatter_maps, material):
    """Run one phonon through the system and record parameters of this run"""

    scattering_types = ScatteringTypes()

    # Run the phonon step-by-step:
    for step_number in range(NUMBER_OF_TIMESTEPS):

        # If phonon has not reached the cold side yet:
        if phonon.is_in_system:

            # Check if different scattering events happened during current time step:
            if INCLUDE_INTERNAL_SCATTERING:
                internal_scattering(phonon, flight, scattering_types)
            surface_scattering(phonon, scattering_types)
            reinitialization(phonon, scattering_types)

            # Record scattering events if any:
            if scattering_types.is_scattered:
                scatter_stats.save_scattering_events(phonon.y, scattering_types)
                flight.add_point_to_path()

            # If no diffuse scattering event occurred, keep measuring the paths and time:
            if not (scattering_types.is_diffuse or scattering_types.is_internal):
                flight.add_step()

            # If diffuse scattering has occurred, reset phonon free path:
            else:
                flight.save_free_paths()
                flight.restart()
                phonon.assign_internal_scattering_time(material)

            # Update scattering and energy maps:
            if OUTPUT_SCATTERING_MAP and scattering_types.is_scattered:
                scatter_maps.add_scattering_to_map(phonon, scattering_types)
            thermal_maps.add_energy_to_maps(phonon, step_number, material)

            # Record time spent in the segment:
            segment_stats.record_time_in_segment(phonon.y)

            # Phonon makes a step forward:
            phonon.move()

            # Reset scattering types for the next step:
            scattering_types.reset()

        # If the phonon reached cold side, record a few parameters and break the loop:
        else:
            flight.finish(step_number)
            break


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
    if not os.path.exists(OUTPUT_FOLDER_NAME):
        os.makedirs(OUTPUT_FOLDER_NAME)
        os.makedirs(OUTPUT_FOLDER_NAME + '/Data')
    shutil.copy('parameters.py', OUTPUT_FOLDER_NAME)
    os.chdir(OUTPUT_FOLDER_NAME)

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

    sys.stdout.write(f'\rSee the results in "{OUTPUT_FOLDER_NAME}" folder.\n')
    sys.stdout.write("\rThank you for using FreePATHS.\n")


if __name__ == "__main__":
    main()
