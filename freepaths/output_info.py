"""Module that outputs general info and basic scattering statistics"""

import time
import logging
import numpy as np
from scipy.constants import electron_volt

from freepaths.config import cf
from freepaths.materials import get_media_class


def output_general_information(start_time, general_stats, mode):
    """This function outputs the simulation information into the Information.txt file"""
    from freepaths.options import SimulationMode
    print(f'\rThe simulation took about {int((time.time() - start_time)//60)} min. to run.')
    n_total = len(general_stats.travel_times)
    percentage = int(100 * np.count_nonzero(general_stats.travel_times) / n_total) if n_total > 0 else 0

    # In MFP sampling mode, NUMBER_OF_PARTICLES is the count per polarization branch,
    # not the true total (see materials.py/_run_branch: each branch is a separate
    # worker looping range(cf.number_of_particles)), so state both explicitly:
    if mode is SimulationMode.PHONON_MFP_SAMPLING:
        particles_line = f'\nNumber of particles = {cf.number_of_particles} per branch ({n_total} total)'
    else:
        particles_line = f'\nNumber of particles = {cf.number_of_particles}'

    info = [
            f'The simulation finished on {time.strftime("%d %B %Y")}, at {time.strftime("%H:%M")}.',
            f'\nIt took about {int((time.time()-start_time)//60)} min to run.\n',
            particles_line,
            f'\nNumber of timesteps = {cf.number_of_timesteps}',
            f'\nLength of a timestep = {cf.timestep} s',
            f'\nTemperature = {cf.temp} K\n',
            f'\nMaterial: {cf.media}\n',
            f'\nLength = {cf.length * 1e9:.1f} nm',
            f'\nWidth = {cf.width * 1e9:.1f} nm',
            f'\nThickness = {cf.thickness * 1e9:.1f} nm\n',
            f'\nSide wall roughness = {cf.side_wall_roughness * 1e9:.1f} nm',
            f'\nTop roughness = {cf.top_roughness * 1e9:.1f} nm',
            f'\nBottom roughness = {cf.bottom_roughness * 1e9:.1f} nm',
            *([ f'\nHole roughness = {cf.hole_roughness * 1e9:.1f} nm'] if cf.holes else []),
            *([ f'\nInterface roughness = {cf.interface_roughness * 1e9:.1f} nm'] if cf.interfaces else []),
            *([ f'\nGrain size = {cf.grain_size * 1e9:.1f} nm, std = {cf.grain_size_std * 1e9:.1f} nm, roughness = {cf.grain_roughness * 1e9:.1f} nm'] if cf.grain_size else []),
            '\n',
            f'\n{percentage:.0f}% of particles reached the cold side\n'
            ]
    with open("Information.txt", "w+", encoding="utf-8") as file:
        file.writelines(info)


def output_scattering_information(scatter_stats):
    """Calculate and output general statistics on scattering events"""

    # Calculate the percentage of different scattering events:
    total_wall = np.sum(scatter_stats.wall_diffuse) + np.sum(scatter_stats.wall_specular)
    total_topbot = np.sum(scatter_stats.top_diffuse) + np.sum(scatter_stats.top_specular)
    total_hole = np.sum(scatter_stats.hole_diffuse) + np.sum(scatter_stats.hole_specular)
    total_pill = np.sum(scatter_stats.pillar_diffuse) + np.sum(scatter_stats.pillar_specular)
    # total_interf counts ALL interface events (both reflection and transmission),
    # because scattering_types.interfaces is set for both. total_transmission is
    # a subset of total_interf, so must NOT be added separately to avoid double-counting.
    total_interf = np.sum(scatter_stats.interfaces_diffuse) + np.sum(scatter_stats.interfaces_specular)
    total_transmission = np.sum(scatter_stats.interfaces_transmission_specular) + np.sum(scatter_stats.interfaces_transmission_diffuse)
    total_retherm = np.sum(scatter_stats.hot_side)
    total_internal = np.sum(scatter_stats.internal)
    total = total_wall + total_topbot + total_hole + total_pill + total_interf + total_retherm + total_internal

    sc_on_walls = 100*(np.sum(scatter_stats.wall_diffuse) +
                         np.sum(scatter_stats.wall_specular)) / total
    sc_on_walls_diff = 100*np.sum(scatter_stats.wall_diffuse) / total_wall
    sc_on_walls_spec = 100*np.sum(scatter_stats.wall_specular) / total_wall

    sc_on_topbot = 100*(np.sum(scatter_stats.top_diffuse) +
                          np.sum(scatter_stats.top_specular)) / total

    if total_topbot != 0:
        sc_on_topbot_diff = 100*np.sum(scatter_stats.top_diffuse) / total_topbot
        sc_on_topbot_spec = 100*np.sum(scatter_stats.top_specular) / total_topbot
    else:
        sc_on_topbot_diff = 0
        sc_on_topbot_spec = 0

    retherm = 100 * total_retherm / total
    internal = 100 * total_internal / total

    info = [
            f'\n{sc_on_walls:.2f}% - scattering on side walls ',
            f'({sc_on_walls_diff:.2f}% - diffuse, ',
            f'{sc_on_walls_spec:.2f}% - specular)',
            f'\n{sc_on_topbot:.2f}% - scattering on top and bottom walls ',
            f'({sc_on_topbot_diff:.2f}% - diffuse, ',
            f'{sc_on_topbot_spec:.2f}% - specular)',
            f'\n{retherm:.2f}% - rethermalization at the hot side',
            f'\n{internal:.2f}% - internal scattering processes',
            ]

    # If scatterers are present, add their information:
    if cf.holes:
        sc_on_holes = 100*(np.sum(scatter_stats.hole_diffuse) +
                             np.sum(scatter_stats.hole_specular)) / total
        sc_on_holes_diff = 100*np.sum(scatter_stats.hole_diffuse) / total_hole
        sc_on_holes_spec = 100*np.sum(scatter_stats.hole_specular) / total_hole
        info.extend([
                    f'\n{sc_on_holes:.2f}% - scattering on hole walls ',
                    f'({sc_on_holes_diff:.2f}% - diffuse, ',
                    f'{sc_on_holes_spec:.2f}% - specular)']
                    )

    if cf.pillars:
        sc_on_pill = 100*(np.sum(scatter_stats.pillar_diffuse) +
                            np.sum(scatter_stats.pillar_specular)) / total
        sc_on_pill_diff = 100*np.sum(scatter_stats.pillar_diffuse) / total_pill
        sc_on_pill_spec = 100*np.sum(scatter_stats.pillar_specular) / total_pill
        info.extend([
                    f'\n{sc_on_pill:.2f}% - scattering on pillar walls ',
                    f'({sc_on_pill_diff:.2f}% - diffuse, ',
                    f'{sc_on_pill_spec:.2f}% - specular)']
                    )

    if cf.interfaces:
        sc_on_interf = 100 * total_interf / total
        sc_on_interf_diff = 100 * np.sum(scatter_stats.interfaces_diffuse) / total_interf
        sc_on_interf_spec = 100 * np.sum(scatter_stats.interfaces_specular) / total_interf
        trans_on_interf = 100 * total_transmission / total_interf
        refl_on_interf  = 100 * (total_interf - total_transmission) / total_interf
        info.extend([
                    f'\n{sc_on_interf:.2f}% - scattering on interfaces ',
                    f'({sc_on_interf_diff:.2f}% - diffuse, ',
                    f'{sc_on_interf_spec:.2f}% - specular), ',
                    f'of which {trans_on_interf:.2f}% - transmission, ',
                    f'{refl_on_interf:.2f}% - reflection',
                    ])

    # Write the file:
    with open("Information.txt", "a", encoding="utf-8") as file:
        file.writelines(info)

    # Warn if multiple scattering events co-occur in a single timestep:
    total_timesteps_with_scattering = np.sum(scatter_stats.total)
    if total_timesteps_with_scattering > 0:
        multi_event_percentage = 100 * (total - total_timesteps_with_scattering) / total_timesteps_with_scattering
        if multi_event_percentage > 5:
            logging.warning(f"In {multi_event_percentage:.1f}% of scattering timesteps, more than one scattering "
                            f"mechanism fired simultaneously. Consider reducing TIMESTEP.")


def output_parameter_warnings(mode, general_stats):
    """Check if parameters used for this simulation made sense considering the simulation results"""

    from freepaths.options import SimulationMode

    # These checks are only relevant for phonon tracing simulations:
    if mode is SimulationMode.PHONON_TRACING:
        travel_times = np.loadtxt("Data/All travel times.csv", encoding='utf-8')

        total_time = cf.timestep * cf.number_of_virtual_timesteps
        time_of_stabilization = cf.number_of_stabilization_timeframes * total_time / cf.number_of_timeframes
        long_travel_times = travel_times[travel_times > time_of_stabilization]
        percentage = (len(long_travel_times) / len(travel_times)) * 100
        if percentage > 10:
            logging.warning(f"Travel time of {percentage}% of particles was longer than the stabilization period.\n" +
                              "Increase stabilization period as the thermal conductivity might be incorrect.")

    # Check how many particles reached the cold side during simulation.
    # Irrelevant in MFP sampling mode, where free paths are measured along the way
    # and reaching the cold side is not required:
    if mode is not SimulationMode.PHONON_MFP_SAMPLING:
        n_total = len(general_stats.travel_times)
        percentage = int(100 * np.count_nonzero(general_stats.travel_times) / n_total) if n_total > 0 else 0
        if percentage < 95:
            message = f"Only {percentage}% of particles reached the cold side. Increase NUMBER_OF_TIMESTEPS."
            # The timeframe bias only concerns thermal profiles, i.e. phonon tracing mode:
            if mode is SimulationMode.PHONON_TRACING:
                message += ("\nParticles removed mid-flight are missing from later timeframes,\n" +
                            "so temperature and heat flux profiles may be biased.")
            logging.warning(message)

    # Pixel size checks are only relevant for phonon tracing mode:
    if mode is SimulationMode.PHONON_TRACING:
        speeds = np.loadtxt("Data/All group velocities.csv", encoding='utf-8')
        if max(speeds) * cf.timestep > cf.length / cf.number_of_pixels_y:
            logging.warning("Pixels in y direction are smaller than length of one step")
        if max(speeds) * cf.timestep > cf.width / cf.number_of_pixels_x:
            logging.warning("Pixels in x direction are smaller than length of one step")

    # Check electron scattering times against timestep:
    if mode is SimulationMode.ELECTRON:
        scattering_times = np.genfromtxt("Data/Scattering time vs energy.csv", delimiter=',', skip_header=1)[:, 1]
        n_short = int(np.sum(scattering_times < cf.timestep))
        percentage = 100 * n_short / len(scattering_times)
        if percentage > 10:
            logging.warning(f"{percentage:.0f}% of scattering times are shorter than TIMESTEP. Consider reducing it.")


def output_electron_information(electron_computations):
    """Print transport properties at the chosen Fermi level."""
    material = get_media_class(cf.media)(cf.temp, fermi_level=cf.media_fermi_level)
    ef_J = material.fermi_level
    ef_meV = ef_J * 1e3 / electron_volt

    def at_ef(array_2col):
        return np.interp(ef_J, array_2col[:, 0], array_2col[:, 1])

    sigma = at_ef(electron_computations.mc_conductivity)

    print(f"For requested Fermi energy of {ef_meV:.1f} meV, σ = {sigma * 1e-3:.2f} kS/m")
    if electron_computations.mc_zt is None:
        print("Phonon thermal conductivity was not pre-computed, so ZT was not calculated.")
