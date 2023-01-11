"""Module that outputs general info and basic scattering statistics"""

import time
import numpy as np

from parameters import *


def output_general_information(start_time):
    """This function outputs the simulation information into the Information.txt file"""
    exit_angles = np.loadtxt("Data/All exit angles.csv")
    percentage = int(100 * np.count_nonzero(exit_angles) / NUMBER_OF_PHONONS)
    print(f'\r{percentage}% of phonons reached the cold side.')
    print(f'The simulation took about {int((time.time() - start_time)//60)} min. to run.')

    with open("Information.txt", "w+") as f:
        info = (
                f'The simulation finished on {time.strftime("%d %B %Y")}, at {time.strftime("%H:%M")}.',
                f'\nIt took about {int((time.time()-start_time)//60)} min to run.\n',
                f'\nNumber of phonons = {NUMBER_OF_PHONONS}',
                f'\nNumber of timesteps = {NUMBER_OF_TIMESTEPS}',
                f'\nLength of a timestep = {TIMESTEP} s',
                f'\nTemperature = {T} K\n',
                f'\nLength = {LENGTH * 1e9:.1f} nm',
                f'\nWidth = {WIDTH * 1e9:.1f} nm',
                f'\nThickness = {THICKNESS * 1e9:.1f} nm\n',
                f'\nSide wall roughness = {SIDE_WALL_ROUGHNESS * 1e9:.1f} nm',
                f'\nHole roughness = {HOLE_ROUGHNESS * 1e9:.1f} nm',
                f'\nTop roughness = {TOP_ROUGHNESS * 1e9:.1f} nm',
                f'\nBottom roughness = {BOTTOM_ROUGHNESS * 1e9:.1f} nm\n',
                f'\n{percentage:.0f}% of phonons reached the cold side\n'
        )
        f.writelines(info)


def output_scattering_information(scatter_stats):
    """Calculate and output general statistics on scattering events"""

    # Calculate the percentage of different scattering events:
    total = np.sum(scatter_stats.total)
    total_wall = np.sum(scatter_stats.wall_diffuse) + np.sum(scatter_stats.wall_specular)
    total_topbot = np.sum(scatter_stats.top_diffuse) + np.sum(scatter_stats.top_specular)
    total_hole = np.sum(scatter_stats.hole_diffuse) + np.sum(scatter_stats.hole_specular)
    total_pill = np.sum(scatter_stats.pillar_diffuse) + np.sum(scatter_stats.pillar_specular)

    scat_on_walls = 100*(np.sum(scatter_stats.wall_diffuse) +
                         np.sum(scatter_stats.wall_specular)) / total
    scat_on_walls_diff = 100*np.sum(scatter_stats.wall_diffuse) / total_wall
    scat_on_walls_spec = 100*np.sum(scatter_stats.wall_specular) / total_wall

    scat_on_topbot = 100*(np.sum(scatter_stats.top_diffuse) +
                          np.sum(scatter_stats.top_specular)) / total
    scat_on_topbot_diff = 100*np.sum(scatter_stats.top_diffuse) / total_topbot
    scat_on_topbot_spec = 100*np.sum(scatter_stats.top_specular) / total_topbot

    if INCLUDE_HOLES:
        scat_on_holes = 100*(np.sum(scatter_stats.hole_diffuse) +
                             np.sum(scatter_stats.hole_specular)) / total
        scat_on_holes_diff = 100*np.sum(scatter_stats.hole_diffuse) / total_hole
        scat_on_holes_spec = 100*np.sum(scatter_stats.hole_specular) / total_hole

    if INCLUDE_PILLARS:
        scat_on_pill = 100*(np.sum(scatter_stats.pillar_diffuse) +
                            np.sum(scatter_stats.pillar_specular)) / total
        scat_on_pill_diff = 100*np.sum(scatter_stats.pillar_diffuse) / total_pill
        scat_on_pill_spec = 100*np.sum(scatter_stats.pillar_specular) / total_pill

    retherm = 100*np.sum(scatter_stats.hot_side) / total
    internal = 100*np.sum(scatter_stats.internal) / total

    # Create output:
    info1 = (
            f'\n{scat_on_walls:.2f}% - scattering on side walls ',
            f'({scat_on_walls_diff:.2f}% - diffuse, ',
            f'{scat_on_walls_spec:.2f}% - specular)',
            f'\n{scat_on_topbot:.2f}% - scattering on top and bottom walls ',
            f'({scat_on_topbot_diff:.2f}% - diffuse, ',
            f'{scat_on_topbot_spec:.2f}% - specular)',
            f'\n{retherm:.2f}% - rethermalization at the hot side',
            f'\n{internal:.2f}% - internal scattering processes',
    )

    if INCLUDE_HOLES:
        info2 = (
                f'\n{scat_on_holes:.2f}% - scattering on hole walls ',
                f'({scat_on_holes_diff:.2f}% - diffuse, ',
                f'{scat_on_holes_spec:.2f}% - specular)',
        )

    if INCLUDE_PILLARS:
        info3 = (
                f'\n{scat_on_pill:.2f}% - scattering on pillar walls ',
                f'({scat_on_pill_diff:.2f}% - diffuse, ',
                f'{scat_on_pill_spec:.2f}% - specular)'
        )

    # Write info into a text file:
    with open("Information.txt", "a") as f:
        f.writelines(info1)
        if INCLUDE_HOLES:
            f.writelines(info2)
        if INCLUDE_PILLARS:
            f.writelines(info3)