"""Main module of the program"""

import os
import sys
import time
import shutil
import colorama
import multiprocessing
import logging
from colorama import Fore, Style

# Modules:
from freepaths.config import cf
from freepaths.run_phonon import run_phonon
from freepaths.phonon import Phonon
from freepaths.flight import Flight
from freepaths.data import ScatteringData, GeneralData, SegmentData, PathData, TriangleScatteringData
from freepaths.progress import Progress
from freepaths.materials import Si, SiC, Graphite
from freepaths.maps import ScatteringMap, ThermalMaps
from freepaths.output_info import output_general_information, output_scattering_information, output_parameter_warnings
from freepaths.animation import create_animation
from freepaths.output_plots import plot_data


class PhononSimulator:
    """
    This class can simulate a number of phonons and save all their data and then return it all
    It is meant to be used as a worker for multiprocessing
    """

    def __init__(self, worker_id, total_phonons, shared_list, output_trajectories_of):

        # Initialize the material:
        if cf.media == "Si":
            self.material = Si(cf.temp)
        elif cf.media == "SiC":
            self.material = SiC(cf.temp)
        elif cf.media == "Graphite":
            self.material = Graphite(cf.temp)
        else:
            logging.error(f"Material {cf.media} is not supported")
            sys.exit()

        # Save some general information about the process:
        self.worker_id = worker_id
        self.total_phonons = total_phonons
        self.result_queue = shared_list
        self.creation_time = time.time()
        self.output_trajectories_of = output_trajectories_of

        # Initiate data structures:
        self.scatter_stats = ScatteringData()
        self.general_stats = GeneralData()
        self.segment_stats = SegmentData()
        self.path_stats = PathData()
        self.places_stats = TriangleScatteringData()
        self.scatter_maps = ScatteringMap()
        self.thermal_maps = ThermalMaps()

        self.total_thermal_conductivity = 0.0

    def simulate_phonon(self, index):
        # Initiate a phonon and its flight:
        phonon = Phonon(self.material)
        flight = Flight(phonon)

        # Run this phonon through the structure:
        run_phonon(phonon, flight, self.scatter_stats, self.places_stats, self.segment_stats, self.thermal_maps, self.scatter_maps, self.material)

        # Record the properties returned for this phonon:
        self.general_stats.save_phonon_data(phonon)
        self.general_stats.save_flight_data(flight)

        # Record trajectories of the first N phonons:
        if index < self.output_trajectories_of:
            self.path_stats.save_phonon_path(flight)

    def simulate_phonons(self, render_progress=False):
        """Simulate a number of phonons and save data to shared datastructure"""

        # Only one of the workers will display its progress as it is similar over all workers:
        if render_progress:
            progress = Progress()

        # Run simulation for each phonon:
        for index in range(self.total_phonons):
            # render progress
            if render_progress:
                progress.render(index, self.total_phonons)

            self.simulate_phonon(index)

        if render_progress:
            progress.render(index+1, self.total_phonons)

        # Collect relevant data:
        collected_data = {
            'scatter_stats': self.scatter_stats.dump_data(),
            'general_stats': self.general_stats.dump_data(),
            'places_stats': self.places_stats.dump_data(),
            'segment_stats': self.segment_stats.dump_data(),
            'path_stats': self.path_stats.dump_data(),
            'scatter_maps': self.scatter_maps.dump_data(),
            'thermal_maps': self.thermal_maps.dump_data(),
            'execution_time': time.time() - self.creation_time,
        }

        # Put the data into shared list:
        self.result_queue.append(collected_data)


def worker_process(worker_id, total_phonons, shared_list, output_trajectories_of, finished_workers):
    try:
        # Create a phononsimulator and run the simulation:
        simulator = PhononSimulator(worker_id, total_phonons, shared_list, output_trajectories_of)
        simulator.simulate_phonons(render_progress=1 if worker_id == 0 else 0)

        # Declare that the calculation is finished:
        finished_workers.value += 1
    except Exception as e:
        sys.stdout.write(f'\rworker {worker_id} had error {e}\n')


def display_workers_finished(finished_workers):
    """ Print out the number of active workers"""
    while True:
        text_to_display = f'  Processes finished: {finished_workers.value}/{cf.num_workers}'
        sys.stdout.write(text_to_display)
        sys.stdout.write(f'\033[{len(text_to_display)}D') # move cursor back
        sys.stdout.flush()
        if finished_workers.value == cf.num_workers: break
        time.sleep(0.3)


def main(input_file):
    """This is the main function, which works under Debye approximation.
    It should be used to simulate phonon paths at low temperatures"""

    sys.stdout.write(f'Simulation of {Fore.GREEN}{cf.output_folder_name}{Style.RESET_ALL}\n')
    start_time = time.time()

    # Create manager for managing variable acces for multiple workers:
    manager = multiprocessing.Manager()

    # These variables created with the manager can safely be accessed by multiple workers:
    shared_list = manager.list()
    finished_workers = manager.Value('i', 0)

    # Divide all the phonons among the workers:
    workload_per_worker = cf.number_of_phonons // cf.num_workers
    remaining_phonons = cf.number_of_phonons % cf.num_workers

    # Divide number of output trajectories to save among workers:
    output_trajectories_per_worker = cf.output_trajectories_of_first // cf.num_workers
    remaining_output_trajectories = cf.output_trajectories_of_first % cf.num_workers

    # Create and start worker processes:
    sys.stdout.write('Starting the workers...\r')
    sys.stdout.flush()
    processes = []
    for i in range(cf.num_workers):
        worker_phonons = workload_per_worker + (1 if i < remaining_phonons else 0)
        output_trajectory_of = output_trajectories_per_worker + (1 if i < remaining_output_trajectories else 0)
        process = multiprocessing.Process(target=worker_process, args=(i, worker_phonons, shared_list, output_trajectory_of, finished_workers))
        processes.append(process)
        process.start()

    # Start a seperate worker to display the number of workers that finished:
    worker_count_process = multiprocessing.Process(target=display_workers_finished, args=(finished_workers,))
    worker_count_process.start()

    # Wait for all processes to finish:
    # Note that join is not called on worker_count_process because we do not want to wait for it to finish
    for process in processes:
        process.join()

    # Wait for the worker count to finish but continue after 3 seconds:
    worker_count_process.join(timeout=3)
    worker_count_process.terminate() # should not be necessary but sometimes process does not terminate

    # Initiate data structures to collect the data from the workers:
    # material = Material(cf.media, num_points=cf.number_of_phonons+1)
    scatter_stats = ScatteringData()
    places_stats = TriangleScatteringData()
    general_stats = GeneralData()
    segment_stats = SegmentData()
    path_stats = PathData()
    scatter_maps = ScatteringMap()
    thermal_maps = ThermalMaps()

    # Collect the results:
    sys.stdout.write('\nCollecting data from workers...\r')

    # Convert the shared list to a normal list so it's easier to use
    result_list = list(shared_list)

    # Check that all workers actually returned some data
    if len(result_list) != cf.num_workers:
        sys.stdout.write(f'WARNING: of {cf.num_workers} workers only the results of {len(result_list)} were collected\n')

    # Put the data from every worker into it's respective place:
    execution_time_list = []
    for collected_data in result_list:
        scatter_stats.read_data(collected_data['scatter_stats'])
        places_stats.read_data(collected_data['places_stats'])
        general_stats.read_data(collected_data['general_stats'])
        segment_stats.read_data(collected_data['segment_stats'])
        path_stats.read_data(collected_data['path_stats'])
        scatter_maps.read_data(collected_data['scatter_maps'])
        thermal_maps.read_data(collected_data['thermal_maps'])
        execution_time_list.append(collected_data['execution_time'])

    # Give some info about the variability in the worker calculation time:
    if cf.num_workers > 1:
        sys.stdout.write(f'Shortest process execution time: {round(min(execution_time_list))}s\n')
        sys.stdout.write(f'Longest process execution time: {round(max(execution_time_list))}s\n')

    # Check if the total number of returned phonons from the workers corresponds with the number of phonons to be simulated:
    if len(general_stats.initial_angles) != cf.number_of_phonons:
        sys.stdout.write(f'WARNING: {cf.number_of_phonons} were meant to be simulated but only {len(general_stats.initial_angles)} phonons were collected from the workers\n')

    # Run additional calculations:
    thermal_maps.calculate_thermal_conductivity()
    thermal_maps.calculate_weighted_flux()
    thermal_maps.calculate_heat_flux_modulus()

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
    sys.stdout.write("\rSaving raw data...")
    general_stats.write_into_files()
    scatter_stats.write_into_files()
    places_stats.write_into_files()
    segment_stats.write_into_files()
    thermal_maps.write_into_files()
    scatter_maps.write_into_files()
    path_stats.write_into_files()

    # Generate animation of phonon paths:
    if cf.output_path_animation:
        create_animation()

    # Analyze and plot the data:
    sys.stdout.write("\rAnalyzing the data...")
    plot_data()

    # Output general information:
    output_general_information(start_time)
    output_scattering_information(scatter_stats)
    output_parameter_warnings()

    sys.stdout.write(f'\rSee the results in {Fore.GREEN}Results/{cf.output_folder_name}{Style.RESET_ALL}\n')
    sys.stdout.write(f"\r{Fore.BLUE}Thank you for using FreePATHS{Style.RESET_ALL}\n\n")
