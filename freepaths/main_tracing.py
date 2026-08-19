"""Main module of the program"""

import os
import sys
import time
import shutil
import signal
import multiprocessing
import traceback
import logging
import numpy as np
from colorama import Fore, Style

# Modules:
from freepaths.config import cf
from freepaths.run_particle import run_particle
from freepaths.electron import Electron, depletion_width
from freepaths.phonon import Phonon
from freepaths.flight import Flight
from freepaths.options import SimulationMode
from freepaths.data import ScatteringData, GeneralData, SegmentData, PathData, TriangleScatteringData
from freepaths.post_computations import ElectronPostComputation
from freepaths.progress import Progress
from freepaths.materials import Si, SiC, Graphite, SiGe
from freepaths.maps import ScatteringMap, ThermalMaps, DriftField
from freepaths.output_info import output_general_information, output_scattering_information, output_parameter_warnings, output_electron_information
from freepaths.animation import create_animation
from freepaths.output_plots import plot_data
from freepaths.scattering import surface_scattering


class ParticleSimulator:
    """
    This class can simulate a number of particles and save all their data and then return it all
    It is meant to be used as a worker for multiprocessing
    """

    def __init__(self, worker_id, mode: SimulationMode, total_particles, shared_list, output_trajectories_of, drift_field=None):

        # Initialize the material:
        if cf.media == "Si":
            self.material = Si(cf.temp)
        elif cf.media == "SiGe":
            self.material = SiGe(cf.temp)
        elif cf.media == "SiC":
            self.material = SiC(cf.temp)
        elif cf.media == "Graphite":
            self.material = Graphite(cf.temp)
        else:
            logging.error(f"Material {cf.media} is not supported")
            sys.exit()

        # Frozen hydrodynamic drift field for this pass (None in the bootstrap pass and
        # in non-hydrodynamic runs); read by momentum-conserving Normal scattering events:
        self.material.drift_field = drift_field

        # Carrier surface-depletion dead-layer width (0 unless SURFACE_POTENTIAL and doping are set):
        cf.depletion_width = depletion_width(self.material)

        # Save some general information about the process:
        self.worker_id = worker_id
        self.total_particles = total_particles
        self.mode = mode
        self.result_queue = shared_list
        self.creation_time = time.time()
        self.output_trajectories_of = output_trajectories_of
        self.number_of_electron_energy_levels = int((cf.energy_upper_bound-cf.energy_lower_bound)/cf.energy_step)

        # Initiate data structures:
        self.scatter_stats = ScatteringData()
        self.general_stats = GeneralData()
        self.segment_stats = SegmentData()
        self.path_stats = PathData()
        self.places_stats = TriangleScatteringData()
        self.scatter_maps = ScatteringMap()
        self.thermal_maps = ThermalMaps() if mode is SimulationMode.PHONON_TRACING else None

        self.total_thermal_conductivity = 0.0
        self.interfaces_transmission_specular = 0
        self.interfaces_transmission_diffuse = 0


    def simulate_particle(self, index):
        # Initiate a particle and its flight:
        if self.mode is SimulationMode.ELECTRON:
            particle = Electron(self.material)
        else:
            particle = Phonon(self.material, self.mode)
        flight = Flight(particle)
        particle.flight = flight

        # Run this particle through the structure:
        run_particle(particle, flight, self.scatter_stats, self.places_stats, self.segment_stats, self.thermal_maps, self.scatter_maps, self.material, self.mode)

        # Record the properties returned for this particle:
        self.general_stats.save_particle_data(particle, self.mode)
        self.general_stats.save_flight_data(flight, self.mode)

        # Record trajectories of the first N particles:
        if index < self.output_trajectories_of:
            self.path_stats.save_particle_path(flight) # FIXME: change method's name

    def simulate_particles(self, render_progress=False):
        """Simulate a number of particles and save data to shared datastructure"""

        # Only one of the workers will display its progress as it is similar over all workers:
        if render_progress:
            progress = Progress()

        # Run simulation for each particle:
        for index in range(self.total_particles):
            # render progress
            if render_progress:
                progress.render(index, self.total_particles)

            self.simulate_particle(index)

        if render_progress:
            progress.render(index+1, self.total_particles)

        # Collect relevant data:
        collected_data = {
            'scatter_stats': self.scatter_stats.dump_data(),
            'general_stats': self.general_stats.dump_data(),
            'places_stats': self.places_stats.dump_data(),
            'segment_stats': self.segment_stats.dump_data(),
            'path_stats': self.path_stats.dump_data(),
            'scatter_maps': self.scatter_maps.dump_data(),
            'thermal_maps': self.thermal_maps.dump_data() if self.thermal_maps is not None else {},
            'execution_time': time.time() - self.creation_time,
        }

        # Put the data into shared list:
        self.result_queue.append(collected_data)


def worker_process(worker_id, mode: SimulationMode, total_particles, shared_list, output_trajectories_of, finished_workers, drift_field=None):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    try:
        # Create a particle simulator and run the simulation:
        simulator = ParticleSimulator(worker_id, mode, total_particles, shared_list, output_trajectories_of, drift_field)
        simulator.simulate_particles(render_progress=1 if worker_id == 0 else 0)

        # Declare that the calculation is finished:
        finished_workers.value += 1
    except Exception as e:
        logging.exception("Worker %d crashed", worker_id)
        # Remonter l’info côté parent
        shared_list.append({
            "worker_id": worker_id,
            "error": traceback.format_exc(),
        })
        # Optionnel : re‑lever pour un exitcode ≠ 0
        sys.exit(1)


def display_workers_finished(finished_workers):
    """ Print out the number of active workers"""
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    try:
        while True:
            text_to_display = f'  Processes finished: {finished_workers.value}/{cf.num_workers}'
            sys.stdout.write(text_to_display)
            sys.stdout.write(f'\033[{len(text_to_display)}D') # move cursor back
            sys.stdout.flush()
            if finished_workers.value == cf.num_workers: break
            time.sleep(0.3)
    except Exception:
        pass


def run_pass(mode: SimulationMode, drift_field=None, total_particles=None):
    """
    Launch all worker processes for one full sweep of the ensemble (one "pass"), with
    the given frozen drift field, and return the collected result_list. In a
    non-hydrodynamic run this is called exactly once; in hydrodynamic mode it is called
    once per self-consistency pass, each time with the drift field frozen at the value
    converged from the previous pass. total_particles overrides cf.number_of_particles
    for this pass (used to run cheaper preruns than the reported run).
    """
    if total_particles is None:
        total_particles = cf.number_of_particles
    manager = multiprocessing.Manager()
    shared_list = manager.list()
    finished_workers = manager.Value('i', 0)

    # Divide all the particles among the workers:
    workload_per_worker = total_particles // cf.num_workers
    remaining_particles = total_particles % cf.num_workers

    # Divide number of output trajectories to save among workers:
    output_trajectories_per_worker = cf.output_trajectories_of_first // cf.num_workers
    remaining_output_trajectories = cf.output_trajectories_of_first % cf.num_workers

    # Create and start worker processes:
    sys.stdout.write('Starting the workers...\r')
    sys.stdout.flush()
    processes = []
    for i in range(cf.num_workers):
        worker_particles = workload_per_worker + (1 if i < remaining_particles else 0)
        output_trajectory_of = output_trajectories_per_worker + (1 if i < remaining_output_trajectories else 0)
        process = multiprocessing.Process(target=worker_process, args=(i, mode, worker_particles, shared_list, output_trajectory_of, finished_workers, drift_field))
        processes.append(process)
        process.start()

    # Start a seperate worker to display the number of workers that finished:
    worker_count_process = multiprocessing.Process(target=display_workers_finished, args=(finished_workers,))
    worker_count_process.start()

    # Wait for all processes to finish:
    # Note that join is not called on worker_count_process because we do not want to wait for it to finish
    try:
        for process in processes:
            process.join()
    except KeyboardInterrupt:
        for process in processes:
            process.terminate()
        raise

    # Wait for the worker count to finish but continue after 3 seconds:
    worker_count_process.join(timeout=3)
    worker_count_process.terminate() # should not be necessary but sometimes process does not terminate
    return list(shared_list)


def merge_results(result_list, mode: SimulationMode):
    """Merge the per-worker data dictionaries from one pass into unified data structures."""
    scatter_stats = ScatteringData()
    places_stats = TriangleScatteringData()
    general_stats = GeneralData()
    segment_stats = SegmentData()
    path_stats = PathData()
    scatter_maps = ScatteringMap()
    thermal_maps = ThermalMaps() if mode is SimulationMode.PHONON_TRACING else None

    if len(result_list) != cf.num_workers:
        sys.stdout.write(f'WARNING: of {cf.num_workers} workers only the results of {len(result_list)} were collected\n')

    execution_time_list = []
    for collected_data in result_list:
        scatter_stats.read_data(collected_data['scatter_stats'])
        places_stats.read_data(collected_data['places_stats'])
        general_stats.read_data(collected_data['general_stats'])
        segment_stats.read_data(collected_data['segment_stats'])
        path_stats.read_data(collected_data['path_stats'])
        scatter_maps.read_data(collected_data['scatter_maps'])
        if thermal_maps is not None:
            thermal_maps.read_data(collected_data['thermal_maps'])
        execution_time_list.append(collected_data['execution_time'])

    return scatter_stats, places_stats, general_stats, segment_stats, path_stats, scatter_maps, thermal_maps, execution_time_list


def main(input_file, mode: SimulationMode):
    """
    This is the main function.
    It should be used to simulate electron paths or phonon paths at low temperatures.
    """

    sys.stdout.write(f'Simulation of {Fore.GREEN}{cf.output_folder_name}{Style.RESET_ALL}\n')
    start_time = time.time()

    # Hydrodynamic mode first traces the ensemble in NUMBER_OF_HYDRODYNAMIC_PRERUNS
    # preruns that only build the self-consistent drift field (frozen within each prerun,
    # under-relaxed between them), then a final reported run that reads the converged field
    # and supplies the reported statistics and maps. The reported run does NOT update the
    # field (there is no pass after it), so unlike a plain "report the last of N passes"
    # scheme the reported run reads the fully converged field. Everything non-hydrodynamic
    # runs a single pass.
    hydrodynamic = cf.phonon_hydrodynamic and mode is SimulationMode.PHONON_TRACING
    number_of_preruns = cf.number_of_hydrodynamic_preruns if hydrodynamic else 0
    number_of_passes = number_of_preruns + 1

    # Material instance for the momentum-susceptibility constant (drift-field derivation):
    hydro_material = {"Si": Si, "SiGe": SiGe, "SiC": SiC, "Graphite": Graphite}[cf.media](cf.temp) if hydrodynamic else None

    drift_field = None
    drift_convergence = []   # per-prerun (n, mean|u_fresh|, mean|u_field|, rel_change) for the convergence test
    flux_profiles = []       # per-pass mid-length cross-width |J_y|(x), to watch it become parabolic
    # Preruns only need to converge the field, so they can use fewer particles than the
    # reported run (NUMBER_OF_HYDRODYNAMIC_PRERUN_PARTICLES; None = same as the reported run):
    prerun_particles = cf.number_of_hydrodynamic_prerun_particles
    for pass_index in range(number_of_passes):
        is_reported_run = (pass_index == number_of_passes - 1)
        pass_particles = None
        if hydrodynamic and not is_reported_run and prerun_particles:
            pass_particles = prerun_particles
        if hydrodynamic:
            n = pass_particles if pass_particles else cf.number_of_particles
            label = "reported run" if is_reported_run else f"prerun {pass_index + 1}/{number_of_preruns}"
            sys.stdout.write(f'\nHydrodynamic {label} ({n} particles)\n')

        result_list = run_pass(mode, drift_field, pass_particles)

        # Collect the results of this pass:
        sys.stdout.write('\nCollecting data from workers...\r')
        (scatter_stats, places_stats, general_stats, segment_stats, path_stats,
         scatter_maps, thermal_maps, execution_time_list) = merge_results(result_list, mode)

        # During preruns, update the self-consistent drift field for the next pass. The
        # freshly measured field is under-relaxed against the previous one (Robbins-Monro)
        # for stability. The reported run still measures its drift (for the saved profile)
        # but does not feed it back:
        if hydrodynamic and thermal_maps is not None:
            thermal_maps.calculate_drift_velocity(hydro_material)

            # Cross-width heat-flux profile at mid-length, one per pass, so the profile can
            # be watched evolving from blunt (bootstrap, no drift) toward parabolic (Poiseuille):
            ny = thermal_maps.heat_flux_map_y.shape[0]
            band = thermal_maps.heat_flux_map_y[int(0.35 * ny):int(0.65 * ny), :]
            flux_profiles.append(np.abs(band).mean(axis=0))

            if not is_reported_run:
                fresh = DriftField(thermal_maps.drift_velocity_x, thermal_maps.drift_velocity_y, thermal_maps.drift_velocity_z)
                previous = drift_field
                if drift_field is None:
                    drift_field = fresh
                else:
                    g = cf.hydrodynamic_preruns_weight
                    drift_field = DriftField(
                        (1 - g) * drift_field.u_x + g * fresh.u_x,
                        (1 - g) * drift_field.u_y + g * fresh.u_y,
                        (1 - g) * drift_field.u_z + g * fresh.u_z,
                    )

                # Convergence diagnostics: mean drift magnitude of the freshly measured
                # field and of the (under-relaxed) field carried to the next prerun, plus
                # the relative L1 change of that field vs the previous prerun. When the
                # relative change flattens toward its noise floor, the preruns are adequate.
                def _mean_mag(f):
                    mag = np.hypot(f.u_x, f.u_y)
                    return float(mag[mag > 0].mean()) if np.any(mag > 0) else 0.0
                if previous is None:
                    rel_change = np.nan
                else:
                    delta = np.abs(drift_field.u_x - previous.u_x) + np.abs(drift_field.u_y - previous.u_y)
                    scale = np.abs(drift_field.u_x) + np.abs(drift_field.u_y)
                    rel_change = float(delta.sum() / (scale.sum() + 1e-30))
                drift_convergence.append((pass_index + 1, _mean_mag(fresh), _mean_mag(drift_field), rel_change))

    # Give some info about the variability in the worker calculation time:
    if cf.num_workers > 1:
        sys.stdout.write(f'Shortest process execution time: {round(min(execution_time_list))}s\n')
        sys.stdout.write(f'Longest process execution time: {round(max(execution_time_list))}s\n')

    # Run thermal calculations on the final pass:
    if thermal_maps is not None:
        thermal_maps.calculate_thermal_conductivity()
        thermal_maps.calculate_heat_flux_modulus()

    # Run calculations on electrons if needed:
    if mode is SimulationMode.ELECTRON:

        electron_computations = ElectronPostComputation(general_stats)
        electron_computations.compute()

    # Create the folder if it does not exist and copy input file there:
    if not os.path.exists(f"Results/{cf.output_folder_name}"):
        os.makedirs(f"Results/{cf.output_folder_name}")
        os.makedirs(f"Results/{cf.output_folder_name}/Data")
    if cf.output_path_animation and not os.path.exists(f"Results/{cf.output_folder_name}/Frames"):
        os.makedirs(f"Results/{cf.output_folder_name}/Frames")
    if input_file:
        shutil.copy(input_file, "Results/" + cf.output_folder_name)
    os.chdir("Results/" + cf.output_folder_name)

    # Save the hydrodynamic-drift-field convergence trace (one row per prerun):
    if drift_convergence:
        np.savetxt("Data/Drift convergence.csv", np.array(drift_convergence), fmt="%1.4e", delimiter=",",
                   header="prerun, mean|u_fresh| (m/s), mean|u_field| (m/s), rel_change_of_field", encoding="utf-8")
    if flux_profiles:
        prof = np.array(flux_profiles).T                       # rows = x pixels, cols = passes
        nx = prof.shape[0]
        x = (np.arange(nx) + 0.5) * cf.width / nx - cf.width / 2
        header = "x (m), " + ", ".join([f"pass{i + 1}" for i in range(prof.shape[1])])
        np.savetxt("Data/Drift convergence profiles.csv", np.column_stack([x, prof]), fmt="%1.4e",
                   delimiter=",", header=header, encoding="utf-8")

    # Try to load pre-computed phonon thermal conductivity to enable ZT calculation:
    if mode is SimulationMode.ELECTRON:
        kappa_ph = None
        if os.path.isfile("Data/Thermal conductivity from MFP.csv"):
            try:
                # The MC electrical conductivity is a solid-convention quantity: the TDF
                # uses the bulk DOS, so holes enter only via the time-of-flight slowdown,
                # which includes both hole-boundary scattering and tortuosity (detours),
                # but not the (1 - porosity) carrier removal. The consistent phonon
                # counterpart is kappa_solid = K_material / (1 + porosity/2), i.e.
                # K_material with the tortuosity part of the Eucken factor applied:
                data = np.atleast_1d(np.genfromtxt("Data/Thermal conductivity from MFP.csv", delimiter=',', skip_header=1)).flatten()
                kappa_ph = float(data[0]) / (1 + float(data[1]) / 2)
            except Exception:
                pass
        elif os.path.isfile("Data/Average thermal conductivity.csv"):
            try:
                data = np.genfromtxt("Data/Average thermal conductivity.csv", delimiter=',', skip_header=1)
                kappa_ph = float(np.atleast_1d(data)[1])  # column 1 = av_mat_tc
            except Exception:
                pass
        if kappa_ph is not None:
            electron_computations.compute_zt(kappa_ph)

    # Save data into files:
    sys.stdout.write("\rSaving raw data...")
    general_stats.write_into_files(mode)
    scatter_stats.write_into_files()
    if cf.holes:
        places_stats.write_into_files()
    segment_stats.write_into_files()
    if thermal_maps is not None:
        thermal_maps.write_into_files(mode)
    if cf.output_scattering_map:
        scatter_maps.write_into_files()
    path_stats.write_into_files()
    if mode is SimulationMode.ELECTRON:
        electron_computations.write_into_file()

    # Generate animation of particle paths:
    if cf.output_path_animation:
        create_animation()

    # Analyze and plot the data:
    sys.stdout.write("\rAnalyzing the data...")
    plot_data(mode, cf)

    # Output general information:
    output_general_information(start_time, general_stats, mode)
    output_scattering_information(scatter_stats)
    if mode is SimulationMode.ELECTRON:
        output_electron_information(electron_computations)
    output_parameter_warnings(mode, general_stats)

    sys.stdout.write(f'\rSee the results in {Fore.GREEN}Results/{cf.output_folder_name}{Style.RESET_ALL}\n')
    sys.stdout.write(f"\r{Fore.BLUE}Thank you for using FreePATHS{Style.RESET_ALL}\n\n")
