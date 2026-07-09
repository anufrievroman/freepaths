"""Module to run thermal conductivity calculations by integrating the phonon dispersion"""

import os
import numpy as np
import sys
import time
import shutil
import scipy
import math
import logging
import signal
import traceback
import multiprocessing
from colorama import Fore, Style

# Modules:
from freepaths.animation import create_animation
from freepaths.config import cf
from freepaths.run_particle import run_particle
from freepaths.phonon import Phonon
from freepaths.flight import Flight
from freepaths.particle_types import ParticleType
from freepaths.data import ScatteringData, GeneralData, SegmentData, PathData, TriangleScatteringData
from freepaths.materials import Si, SiC, Graphite, SiGe
from freepaths.maps import ScatteringMap
from freepaths.output_info import output_general_information, output_scattering_information, output_parameter_warnings
from freepaths.output_plots import plot_data


def _branch_worker(branch_number, shared_list):
    """Process all phonons for one polarization branch and return results via shared_list."""
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    try:
        _run_branch(branch_number, shared_list)
    except Exception:
        logging.exception("Branch worker %d crashed", branch_number)


def _run_branch(branch_number, shared_list):

    if cf.media == "Si":
        material = Si(cf.temp, num_points=cf.number_of_particles + 1)
    elif cf.media == "SiGe":
        material = SiGe(cf.temp, num_points=cf.number_of_particles + 1)
    elif cf.media == "SiC":
        material = SiC(cf.temp, num_points=cf.number_of_particles + 1)
    elif cf.media == "Graphite":
        material = Graphite(cf.temp, num_points=cf.number_of_particles + 1)
    else:
        logging.error(f"Material {cf.media} is not supported")
        return

    scatter_stats = ScatteringData()
    general_stats = GeneralData()
    segment_stats = SegmentData()
    places_stats = TriangleScatteringData()
    path_stats = PathData()
    scatter_maps = ScatteringMap()

    total_thermal_conductivity = 0.0

    for index in range(cf.number_of_particles):
        k_vector = (material.dispersion[index + 1, 0] + material.dispersion[index, 0]) / 2
        d_k_vector = material.dispersion[index + 1, 0] - material.dispersion[index, 0]

        phonon = Phonon(material, branch_number, index)
        flight = Flight(phonon)

        run_particle(phonon, flight, scatter_stats, places_stats, segment_stats, None, scatter_maps, material)

        omega = 2 * math.pi * phonon.f
        part = scipy.constants.hbar * omega / (scipy.constants.k * cf.temp)
        c_p = scipy.constants.k * part**2 * math.exp(part) / (math.exp(part) - 1)**2

        mean_relax_time = flight.mean_free_path / phonon.speed
        flight.thermal_conductivity = (1 / (6 * math.pi**2)) * c_p * phonon.speed**2 * mean_relax_time * k_vector**2 * d_k_vector
        total_thermal_conductivity += flight.thermal_conductivity

        general_stats.save_particle_data(phonon)
        general_stats.save_flight_data(flight)

        if index < cf.output_trajectories_of_first:
            path_stats.save_particle_path(flight)

    shared_list.append({
        'total_thermal_conductivity': total_thermal_conductivity,
        'scatter_stats': scatter_stats.dump_data(),
        'general_stats': general_stats.dump_data(),
        'segment_stats': segment_stats.dump_data(),
        'places_stats': places_stats.dump_data(),
        'path_stats': path_stats.dump_data(),
        'scatter_maps': scatter_maps.dump_data(),
    })


def main(input_file, particle_type):
    """This is the main function, which integrates phonon dispersion to get thermal conductivity"""

    print(f'Mean free path sampling of {Fore.GREEN}{cf.output_folder_name}{Style.RESET_ALL}')
    start_time = time.time()

    # Spawn one worker process per polarization branch:
    manager = multiprocessing.Manager()
    shared_list = manager.list()

    processes = [
        multiprocessing.Process(target=_branch_worker, args=(branch_number, shared_list))
        for branch_number in range(3)
    ]
    for process in processes:
        process.start()
    try:
        for process in processes:
            process.join()
    except KeyboardInterrupt:
        for process in processes:
            process.terminate()
        raise

    # Collect and merge results from all three branches:
    scatter_stats = ScatteringData()
    general_stats = GeneralData()
    segment_stats = SegmentData()
    places_stats = TriangleScatteringData()
    path_stats = PathData()
    scatter_maps = ScatteringMap()
    total_thermal_conductivity = 0.0

    for result in shared_list:
        total_thermal_conductivity += result['total_thermal_conductivity']
        scatter_stats.read_data(result['scatter_stats'])
        general_stats.read_data(result['general_stats'])
        segment_stats.read_data(result['segment_stats'])
        places_stats.read_data(result['places_stats'])
        path_stats.read_data(result['path_stats'])
        scatter_maps.read_data(result['scatter_maps'])

    # Create the output folder and copy input file:
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
    if cf.output_scattering_map:
        scatter_maps.write_into_files()
    path_stats.write_into_files()

    if cf.output_path_animation:
        create_animation()

    np.savetxt("Data/Thermal conductivity from MFP.csv", np.array([total_thermal_conductivity]), fmt='%2.4e', header="K [W/mK]", encoding='utf-8')

    sys.stdout.write("\rAnalyzing the data...")
    plot_data(particle_type, cf, mfp_sampling=True)

    output_general_information(start_time)
    output_scattering_information(scatter_stats)
    output_parameter_warnings(ParticleType.PHONON)
    sys.stdout.write(f'\rSee the results in {Fore.GREEN}Results/{cf.output_folder_name}{Style.RESET_ALL}\n')
    sys.stdout.write(f"\rThermal conductivity = {Fore.GREEN}{total_thermal_conductivity:.5f}{Style.RESET_ALL} W/m·K\n")
    sys.stdout.write(f"\r{Fore.BLUE}Thank you for using FreePATHS{Style.RESET_ALL}\n\n")
