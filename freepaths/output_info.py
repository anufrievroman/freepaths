"""Module that outputs general info and basic scattering statistics"""

import time
import logging
import numpy as np

from colorama import Fore, Style

from freepaths.config import cf


def output_general_information(start_time):
    """This function outputs the simulation information into the Information.txt file"""
    print(f'\rThe simulation took about {int((time.time() - start_time)//60)} min. to run.')
    travel_times = np.loadtxt("Data/All travel times.csv", encoding='utf-8')
    percentage = int(100 * np.count_nonzero(travel_times) / cf.number_of_phonons)

    info = [
            f'The simulation finished on {time.strftime("%d %B %Y")}, at {time.strftime("%H:%M")}.',
            f'\nIt took about {int((time.time()-start_time)//60)} min to run.\n',
            f'\nNumber of phonons = {cf.number_of_phonons}',
            f'\nNumber of timesteps = {cf.number_of_timesteps}',
            f'\nLength of a timestep = {cf.timestep} s',
            f'\nTemperature = {cf.temp} K\n',
            f'\nMaterial: {cf.media}\n',
            f'\nLength = {cf.length * 1e9:.1f} nm',
            f'\nWidth = {cf.width * 1e9:.1f} nm',
            f'\nThickness = {cf.thickness * 1e9:.1f} nm\n',
            f'\nSide wall roughness = {cf.side_wall_roughness * 1e9:.1f} nm',
            f'\nHole roughness = {cf.hole_roughness * 1e9:.1f} nm',
            f'\nTop roughness = {cf.top_roughness * 1e9:.1f} nm',
            f'\nBottom roughness = {cf.bottom_roughness * 1e9:.1f} nm',
            f'\nInterface roughness = {cf.interface_roughness * 1e9:.1f} nm\n',
            f'\n{percentage:.0f}% of phonons reached the cold side\n'
            ]
    with open("Information.txt", "w+", encoding="utf-8") as file:
        file.writelines(info)


def output_scattering_information(scatter_stats):
    """Calculate and output general statistics on scattering events"""

    # Calculate the percentage of different scattering events:
    total = np.sum(scatter_stats.total)
    total_wall = np.sum(scatter_stats.wall_diffuse) + np.sum(scatter_stats.wall_specular)
    total_topbot = np.sum(scatter_stats.top_diffuse) + np.sum(scatter_stats.top_specular)
    total_hole = np.sum(scatter_stats.hole_diffuse) + np.sum(scatter_stats.hole_specular)
    total_pill = np.sum(scatter_stats.pillar_diffuse) + np.sum(scatter_stats.pillar_specular)
    total_interf = np.sum(scatter_stats.interfaces_diffuse) + np.sum(scatter_stats.interfaces_specular)

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

    retherm = 100*np.sum(scatter_stats.hot_side) / total
    internal = 100*np.sum(scatter_stats.internal) / total

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

    if cf. interfaces:
        sc_on_interf = 100*(np.sum(scatter_stats.interfaces_diffuse) +
                             np.sum(scatter_stats.interfaces_specular)) / total
        sc_on_interf_diff = 100*np.sum(scatter_stats.interfaces_diffuse) / total_interf
        sc_on_interf_spec = 100*np.sum(scatter_stats.interfaces_specular) / total_interf
        info.extend([
                    f'\n{sc_on_interf:.2f}% - scattering on interfaces ',
                    f'({sc_on_interf_diff:.2f}% - diffuse, ',
                    f'{sc_on_interf_spec:.2f}% - specular)']
                    )

    # Write the file:
    with open("Information.txt", "a", encoding="utf-8") as file:
        file.writelines(info)


def output_parameter_warnings():
    """Check if parameters used for this simulation made sense considering the simulation results"""

    # Check if some phonons had longer travel times than stabilization period:
    travel_times = np.loadtxt("Data/All travel times.csv", encoding='utf-8')
    total_time = cf.timestep * cf.number_of_timesteps
    time_of_stabilization = cf.number_of_stabilization_timeframes * total_time / cf. number_of_timeframes
    long_travel_times = travel_times[travel_times > time_of_stabilization]
    percentage = (len(long_travel_times) / len(travel_times)) * 100
    if percentage > 10:
        logging.warning(f"Travel time of {percentage}% of phonons was longer than the stabilization period.\n" +
                          "Increase stabilization period as the thermal conductivity might be incorrect.")

    # Check if pixel size is too small:
    speeds = np.loadtxt("Data/All group velocities.csv", encoding='utf-8')
    if max(speeds) * cf.timestep > cf.length / cf.number_of_pixels_y:
        logging.warning("Pixels in y direction are smaller than length of one step")
    if max(speeds) * cf.timestep > cf.width / cf.number_of_pixels_x:
        logging.warning("Pixels in x direction are smaller than length of one step")

    # Check how many phonons reached the cold side during simulation:
    percentage = int(100 * np.count_nonzero(travel_times) / cf.number_of_phonons)
    if percentage < 95:
        logging.warning(f"Only {percentage}% of phonons reached the cold side. Increase number of timesteps.")
