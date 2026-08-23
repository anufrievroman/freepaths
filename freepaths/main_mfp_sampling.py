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
from freepaths.options import SimulationMode
from freepaths.data import ScatteringData, GeneralData, SegmentData, PathData, TriangleScatteringData
from freepaths.materials import Si, SiC, Graphite, SiGe
from freepaths.maps import ScatteringMap
from freepaths.output_info import output_general_information, output_scattering_information, output_parameter_warnings
from freepaths.materials import get_media_class
from freepaths.output_plots import plot_data


def display_all_progress(shared_progress, finished_branches, branch_names):
    """Render a live progress percentage for each branch on its own line."""
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    n = len(branch_names)
    width = max(len(name) for name in branch_names)
    try:
        for name in branch_names:
            sys.stdout.write(f'  {name.ljust(width)}:   0%\n')
        sys.stdout.flush()
        while True:
            sys.stdout.write(f'\033[{n}A')
            for i, name in enumerate(branch_names):
                sys.stdout.write(f'\r  {name.ljust(width)}: {shared_progress[i]:3d}%\n')
            sys.stdout.flush()
            if finished_branches.value == n:
                break
            time.sleep(0.1)
    except Exception:
        pass


def _branch_worker(branch_number, shared_list, finished_branches, shared_progress):
    """Entry point for each worker process. Catches exceptions so they are logged
    rather than silently swallowed by the multiprocessing machinery."""
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    try:
        _run_branch(branch_number, shared_list, shared_progress)
        shared_progress[branch_number] = 100
        finished_branches.value += 1
    except Exception:
        logging.exception("Branch worker %d crashed", branch_number)


def _run_branch(branch_number, shared_list, shared_progress):
    """Run all phonons for one polarization branch and push results to shared_list.

    Each call processes NUMBER_OF_PARTICLES phonons sampled uniformly across the
    dispersion of the given branch.  The k-vector integration weight (d_k_vector)
    is the width of the k-space interval assigned to each phonon, so the sum of
    flight.thermal_conductivity over all phonons approximates the branch integral
    in Phys. Rev. 132 2461 (1963).

    Branch 2 (TA) drives the shared progress bar; the other two run silently.
    Results are appended to shared_list as a dict of dump_data() payloads so that
    the main process can merge them with the standard read_data() mechanism.
    """

    # Each worker builds its own material instance so the dispersion is
    # sampled with the correct number of points (number_of_particles + 1
    # intervals → number_of_particles midpoints):
    if cf.media == "Si":
        material = Si(cf.temp, num_points=cf.number_of_particles + 1)
    elif cf.media == "SiGe":
        material = SiGe(cf.temp, num_points=cf.number_of_particles + 1)
    elif cf.media == "SiC":
        material = SiC(cf.temp, num_points=cf.number_of_particles + 1)
    elif cf.media == "Graphite":
        material = Graphite(cf.temp, num_points=cf.number_of_particles + 1, isotope_c13_concentration=cf.isotope_c13_concentration)
    else:
        logging.error(f"Material {cf.media} is not supported")
        return

    # Local data structures — merged into the global ones after all workers finish:
    scatter_stats = ScatteringData()
    general_stats = GeneralData()
    segment_stats = SegmentData()
    places_stats = TriangleScatteringData()
    path_stats = PathData()
    scatter_maps = ScatteringMap()

    total_thermal_conductivity = 0.0
    total_kappa_variance = 0.0
    total_kappa1 = 0.0        # Callaway kappa1 accumulator (uses combined tau_C)
    total_k2_num = 0.0        # Callaway kappa2 numerator sum:   w * tau_C / tau_N
    total_k2_den = 0.0        # Callaway kappa2 denominator sum:  w * tau_C / (tau_N * tau_R)
    last_pct = -1

    for index in range(cf.number_of_particles):
        pct = 100 * index // cf.number_of_particles
        if pct > last_pct:
            shared_progress[branch_number] = pct
            last_pct = pct

        # Mid-point k-vector and integration weight for this phonon:
        k_vector = (material.dispersion[index + 1, 0] + material.dispersion[index, 0]) / 2
        d_k_vector = material.dispersion[index + 1, 0] - material.dispersion[index, 0]

        # Initialise phonon at this k-point and run it through the geometry:
        phonon = Phonon(material, SimulationMode.PHONON_MFP_SAMPLING, branch_number, index)
        flight = Flight(phonon)
        run_particle(phonon, flight, scatter_stats, places_stats, segment_stats, None, scatter_maps, material, SimulationMode.PHONON_MFP_SAMPLING)

        # Mode heat capacity (Ref. PRB 88 155318, 2013):
        omega = 2 * math.pi * phonon.f
        part = scipy.constants.hbar * omega / (scipy.constants.k * cf.temp)
        c_p = scipy.constants.k * part**2 * math.exp(part) / (math.exp(part) - 1)**2

        # Thermal conductivity contribution from this k-point (Ref. Phys. Rev. 132 2461, 1963).
        # The measured mean free path is limited by resistive (Umklapp + isotope) and boundary
        # scattering, so mean_relax_time is the Callaway resistive time tau_R (Normal scattering
        # does not fire in the MFP-sampling mode). weight is the common per-mode integrand factor.
        mean_relax_time = flight.mean_free_path / phonon.speed
        weight = (1 / (6 * math.pi**2)) * c_p * phonon.speed**2 * k_vector**2 * d_k_vector
        flight.thermal_conductivity = weight * mean_relax_time          # single-mode RTA (SMRT)
        total_thermal_conductivity += flight.thermal_conductivity

        # Callaway dual-relaxation correction (Callaway, Phys. Rev. 113, 1046 (1959)):
        # kappa = kappa1 + kappa2, with 1/tau_C = 1/tau_R + 1/tau_N. kappa1 uses the combined
        # tau_C (Normal provisionally counted as resistive); kappa2 is the momentum-recovery term
        # kappa2 = (sum w*tau_C/tau_N)^2 / (sum w*tau_C/(tau_N*tau_R)), accumulated linearly here
        # and combined in main(). tau_C/tau_N = tau_C*normal_rate. Only computed when
        # PHONON_HYDRODYNAMIC is on (the Normal / momentum-recovery treatment); otherwise
        # normal_rate = 0 and kappa1 collapses to the standard SMRT conductivity.
        normal_rate = material.phonon_normal_rate(omega) if cf.phonon_hydrodynamic else 0.0
        if normal_rate > 0 and mean_relax_time > 0:
            tau_C = 1.0 / (1.0 / mean_relax_time + normal_rate)
            total_kappa1 += weight * tau_C
            total_k2_num += weight * tau_C * normal_rate
            total_k2_den += weight * tau_C * normal_rate / mean_relax_time
        else:
            total_kappa1 += weight * mean_relax_time

        # This mode's kappa is linear in its mean free path, so its variance
        # propagates directly through the same proportionality (kappa_i / mfp_i):
        if flight.mean_free_path > 0:
            kappa_relative_sem = flight.mean_free_path_sem / flight.mean_free_path
            total_kappa_variance += (flight.thermal_conductivity * kappa_relative_sem) ** 2

        general_stats.save_particle_data(phonon, SimulationMode.PHONON_MFP_SAMPLING)
        general_stats.save_flight_data(flight, SimulationMode.PHONON_MFP_SAMPLING)

        if index < cf.output_trajectories_of_first:
            path_stats.save_particle_path(flight)

    # Push this branch's results to the shared list for the main process to collect:
    shared_list.append({
        'total_thermal_conductivity': total_thermal_conductivity,
        'total_kappa_variance': total_kappa_variance,
        'total_kappa1': total_kappa1,
        'total_k2_num': total_k2_num,
        'total_k2_den': total_k2_den,
        'scatter_stats': scatter_stats.dump_data(),
        'general_stats': general_stats.dump_data(),
        'segment_stats': segment_stats.dump_data(),
        'places_stats': places_stats.dump_data(),
        'path_stats': path_stats.dump_data(),
        'scatter_maps': scatter_maps.dump_data(),
    })


def calculate_porosity(grid_points=300):
    """Calculate the porosity (hole area fraction) of the structure by sampling
    the hole shapes on a regular grid. Holes partially outside the simulation
    domain are accounted for with their inside part only."""
    if not cf.holes:
        return 0.0
    xs = np.linspace(-cf.width / 2, cf.width / 2, grid_points)
    ys = np.linspace(0, cf.length, grid_points)
    inside = 0
    for x in xs:
        for y in ys:
            inside += any(hole.is_inside(x, y, None, cf) for hole in cf.holes)
    return inside / (grid_points * grid_points)


def main(input_file, mode: SimulationMode):
    """Integrate the phonon dispersion over all three branches in parallel to get
    the bulk thermal conductivity via the MFP-sampling method."""

    print(f'Mean free path sampling of {Fore.GREEN}{cf.output_folder_name}{Style.RESET_ALL}')
    start_time = time.time()

    branch_names = get_media_class(cf.media).dispersion_branch_names
    n_branches = len(branch_names)

    # Spawn one worker process per polarization branch, all running concurrently:
    manager = multiprocessing.Manager()
    shared_list = manager.list()
    finished_branches = manager.Value('i', 0)
    shared_progress = manager.list([0] * n_branches)

    processes = [
        multiprocessing.Process(target=_branch_worker, args=(branch_number, shared_list, finished_branches, shared_progress))
        for branch_number in range(n_branches - 1, -1, -1)
    ]
    for process in processes:
        process.start()

    display_process = multiprocessing.Process(target=display_all_progress, args=(shared_progress, finished_branches, branch_names))
    display_process.start()

    try:
        for process in processes:
            process.join()
    except KeyboardInterrupt:
        for process in processes:
            process.terminate()
        raise

    display_process.join(timeout=3)
    display_process.terminate()

    # Merge per-branch results into single data structures:
    scatter_stats = ScatteringData()
    general_stats = GeneralData()
    segment_stats = SegmentData()
    places_stats = TriangleScatteringData()
    path_stats = PathData()
    scatter_maps = ScatteringMap()
    total_thermal_conductivity = 0.0
    total_kappa_variance = 0.0
    total_kappa1 = 0.0
    total_k2_num = 0.0
    total_k2_den = 0.0

    for result in shared_list:
        total_thermal_conductivity += result['total_thermal_conductivity']
        total_kappa_variance += result['total_kappa_variance']
        total_kappa1 += result['total_kappa1']
        total_k2_num += result['total_k2_num']
        total_k2_den += result['total_k2_den']
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

    # Save data and generate plots:
    general_stats.write_into_files(mode)
    scatter_stats.write_into_files()
    segment_stats.write_into_files()
    if cf.output_scattering_map:
        scatter_maps.write_into_files()
    path_stats.write_into_files()

    if cf.output_path_animation:
        create_animation()

    # The MFP-sampling integral uses the bulk DoS with boundary-limited free paths,
    # so total_thermal_conductivity is the conductivity of the solid material
    # (MFP suppression only, no volume removed by the holes). To obtain the
    # effective conductivity of the porous structure, it is scaled by the Eucken
    # factor F = (1 - phi) / (1 + phi/2), where phi is the porosity:
    porosity = calculate_porosity()
    eucken_factor = (1 - porosity) / (1 + porosity / 2)
    effective_thermal_conductivity = total_thermal_conductivity * eucken_factor

    header = "K_material [W/mK], Porosity, Eucken factor, K_effective [W/mK]"
    np.savetxt("Data/Thermal conductivity from MFP.csv",
               np.array([[total_thermal_conductivity, porosity, eucken_factor, effective_thermal_conductivity]]),
               fmt='%2.4e', delimiter=',', header=header, encoding='utf-8')

    # Kept in its own file rather than appended to "Thermal conductivity from MFP.csv",
    # since several plotting scripts index into that file positionally (incl. by
    # [-1] for the last column), which a trailing column would silently break:
    kappa_sem = total_kappa_variance ** 0.5
    effective_kappa_sem = kappa_sem * eucken_factor
    np.savetxt("Data/Thermal conductivity from MFP uncertainty.csv",
               np.array([[kappa_sem, effective_kappa_sem]]),
               fmt='%2.4e', delimiter=',',
               header="K_material SEM [W/mK], K_effective SEM [W/mK]", encoding='utf-8')

    # Callaway dual-relaxation conductivity kappa = kappa1 + kappa2 (momentum-recovery term).
    # Only written in hydrodynamic mode (the Normal-process treatment); in the standard RTA the
    # -s output is unchanged. Kept in a separate file so the positional indexers of
    # "Thermal conductivity from MFP.csv" are not disturbed. See main_mfp_sampling._run_branch.
    kappa2 = kappa_callaway = None
    if cf.phonon_hydrodynamic:
        kappa2 = (total_k2_num ** 2 / total_k2_den) if total_k2_den > 0 else 0.0
        kappa_callaway = total_kappa1 + kappa2
        np.savetxt("Data/Thermal conductivity Callaway.csv",
                   np.array([[total_kappa1, kappa2, kappa_callaway]]),
                   fmt='%2.4e', delimiter=',',
                   header="kappa1 [W/mK], kappa2 [W/mK], kappa_Callaway [W/mK]", encoding='utf-8')

    sys.stdout.write("\rAnalyzing the data...")
    plot_data(mode, cf)

    output_general_information(start_time, general_stats, mode)
    output_scattering_information(scatter_stats)
    output_parameter_warnings(mode, general_stats)
    sys.stdout.write(f'\rSee the results in {Fore.GREEN}Results/{cf.output_folder_name}{Style.RESET_ALL}\n')
    if cf.phonon_hydrodynamic:
        sys.stdout.write(f"\rMaterial thermal conductivity (SMRT) = {Fore.GREEN}{total_thermal_conductivity:.5f} ± {kappa_sem:.5f}{Style.RESET_ALL} W/m·K\n")
        sys.stdout.write(f"\rMaterial thermal conductivity (Callaway) = {Fore.GREEN}{kappa_callaway:.5f}{Style.RESET_ALL} W/m·K "
                         f"(kappa1 {total_kappa1:.2f} + kappa2 {kappa2:.2f})\n")
    else:
        sys.stdout.write(f"\rMaterial thermal conductivity = {Fore.GREEN}{total_thermal_conductivity:.5f} ± {kappa_sem:.5f}{Style.RESET_ALL} W/m·K\n")
    if cf.holes:
        sys.stdout.write(f"\rEffective thermal conductivity = {Fore.GREEN}{effective_thermal_conductivity:.5f} ± {effective_kappa_sem:.5f}{Style.RESET_ALL} W/m·K "
                         f"(porosity {porosity:.1%}, Eucken factor {eucken_factor:.3f})\n")
    sys.stdout.write(f"\r{Fore.BLUE}Thank you for using FreePATHS{Style.RESET_ALL}\n\n")
