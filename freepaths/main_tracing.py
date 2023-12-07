"""Main module of the program"""

import os
import sys
import time
import shutil
import colorama
import multiprocessing
from colorama import Fore, Style

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


class PhononSimulator:
    def __init__(self, worker_id, total_phonons, shared_list, output_trajectories_of):
        self.worker_id = worker_id
        self.total_phonons = total_phonons
        self.result_queue = shared_list
        self.creation_time = time.time()
        self.output_trajectories_of = output_trajectories_of
        
        # Initiate data structures:
        self.material = Material(cf.media, num_points=cf.number_of_phonons+1)   # ok
        self.scatter_stats = ScatteringData()                                   # ok
        self.general_stats = GeneralData()                                      # ok
        self.segment_stats = SegmentData()                                      # ok
        self.path_stats = PathData()                                            # ok
        self.scatter_maps = ScatteringMap()                                     # ok
        self.thermal_maps = ThermalMaps()                                       # ok
        
        self.total_thermal_conductivity = 0.0
    
    def simulate_phonons(self, render_progress=False):
        
        if render_progress:
            progress = Progress()
        
        for index in range(self.total_phonons):
            if render_progress:
                progress.render(index, self.total_phonons)
            
            # Initiate a phonon and its flight:
            phonon = Phonon(self.material)
            flight = Flight(phonon)
            
            # Run this phonon through the structure:
            run_phonon(phonon, flight, self.scatter_stats, self.segment_stats, self.thermal_maps, self.scatter_maps, self.material)

            # Record the properties returned for this phonon:
            self.general_stats.save_phonon_data(phonon)
            self.general_stats.save_flight_data(flight)
            
            # Record trajectories of the first N phonons:
            if index < self.output_trajectories_of:
                self.path_stats.save_phonon_path(flight)
        
        # collect relevant data
        collected_data = {
            'scatter_stats': self.scatter_stats.dump_data(),
            'general_stats': self.general_stats.dump_data(),
            'segment_stats': self.segment_stats.dump_data(),
            'path_stats': self.path_stats.dump_data(),
            'scatter_maps': self.scatter_maps.dump_data(),
            'thermal_maps': self.thermal_maps.dump_data(),
            'execution_time': time.time() - self.creation_time,
        }
        
        # put data into queue
        self.result_queue.append(collected_data)


def worker_process(worker_id, total_phonons, shared_list, output_trajectories_of, finished_workers):
    try:
        simulator = PhononSimulator(worker_id, total_phonons, shared_list, output_trajectories_of)
        simulator.simulate_phonons(render_progress=1 if worker_id == 0 else 0)
    except Exception as e:
        print(f'worker {worker_id} had error {e}')
    
    finished_workers.value += 1


def main(input_file):
    """This is the main function, which works under Debye approximation.
    It should be used to simulate phonon paths at low temperatures"""

    print(f'Simulation of {Fore.GREEN}{cf.output_folder_name}{Style.RESET_ALL}')
    start_time = time.time()
    
    # create queue for collection of data
    manager = multiprocessing.Manager()
    shared_list = manager.list()
    
    # Create a shared variable to keep track of workers that are finished
    finished_workers = manager.Value('i', 0)
    
    # Divide workload among workers
    workload_per_worker = cf.number_of_phonons // cf.num_workers
    remaining_phonons = cf.number_of_phonons % cf.num_workers

    output_trajectories_per_worker = cf.output_trajectories_of_first // cf.num_workers
    remaining_output_trajectories = cf.output_trajectories_of_first % cf.num_workers

    # Create and start worker processes
    print('Starting the workers')
    processes = []
    for i in range(cf.num_workers):
        worker_phonons = workload_per_worker + (1 if i < remaining_phonons else 0)
        output_trajectory_of = output_trajectories_per_worker + (1 if i < remaining_output_trajectories else 0)
        process = multiprocessing.Process(target=worker_process, args=(i, worker_phonons, shared_list, output_trajectory_of, finished_workers))
        processes.append(process)
        process.start()
    
    # display number of active workers
    while finished_workers.value != cf.num_workers:
        text_to_display = f'    Workers finished: {finished_workers.value}/{cf.num_workers}'
        sys.stdout.write(text_to_display)
        sys.stdout.write(f'\033[{len(text_to_display)}D') # move cursor back
        sys.stdout.flush()
        time.sleep(0.3)
    
    # Wait for all processes to finish
    for process in processes:
        process.join()
    
    # Initiate data structures:
    # material = Material(cf.media, num_points=cf.number_of_phonons+1)
    scatter_stats = ScatteringData()
    general_stats = GeneralData()
    segment_stats = SegmentData()
    path_stats = PathData()
    scatter_maps = ScatteringMap()
    thermal_maps = ThermalMaps()
    
    # collect the results
    print(f'Collecting data from workers:')
    result_list = list(shared_list)
    while result_list:
        collected_data = result_list.pop()
        scatter_stats.read_data(collected_data['scatter_stats'])
        general_stats.read_data(collected_data['general_stats'])
        segment_stats.read_data(collected_data['segment_stats'])
        path_stats.read_data(collected_data['path_stats'])
        scatter_maps.read_data(collected_data['scatter_maps'])
        thermal_maps.read_data(collected_data['thermal_maps'])
        print(collected_data['execution_time'])

    # Run additional calculations:
    thermal_maps.calculate_thermal_conductivity()
    thermal_maps.calculate_normalized_flux()
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

    # Generate animation of phonon paths:
    if cf.output_path_animation:
        create_animation()

    # Analyze and plot the data:
    sys.stdout.write("\rAnalyzing the data...")
    plot_data()

    # Output general information:
    output_general_information(start_time)
    output_scattering_information(scatter_stats)

    sys.stdout.write(f'\rSee the results in {Fore.GREEN}Results/{cf.output_folder_name}{Style.RESET_ALL}\n')
    sys.stdout.write(f"\r{Fore.BLUE}Thank you for using FreePATHS{Style.RESET_ALL}\n")
