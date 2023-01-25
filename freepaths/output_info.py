"""Module that outputs general info and basic scattering statistics"""

import time
import numpy as np

from freepaths.config import cf


def output_general_information(start_time):
    """This function outputs the simulation information into the Information.txt file"""
    exit_angles = np.loadtxt("Data/All exit angles.csv")
    percentage = int(100 * np.count_nonzero(exit_angles) / cf.number_of_phonons)
    print(f'\r{percentage}% of phonons reached the cold side.')
    print(f'The simulation took about {int((time.time() - start_time)//60)} min. to run.')

    with open("Information.txt", "w+", encoding="utf-8") as file:
        info = (
                f'The simulation finished on {time.strftime("%d %B %Y")}, at {time.strftime("%H:%M")}.',
                f'\nIt took about {int((time.time()-start_time)//60)} min to run.\n',
                f'\nNumber of phonons = {cf.number_of_phonons}',
                f'\nNumber of timesteps = {cf.number_of_timesteps}',
                f'\nLength of a timestep = {cf.timestep} s',
                f'\nTemperature = {cf.temp} K\n',
                f'\nLength = {cf.length * 1e9:.1f} nm',
                f'\nWidth = {cf.width * 1e9:.1f} nm',
                f'\nThickness = {cf.thickness * 1e9:.1f} nm\n',
                f'\nSide wall roughness = {cf.side_wall_roughness * 1e9:.1f} nm',
                f'\nHole roughness = {cf.hole_roughness * 1e9:.1f} nm',
                f'\nTop roughness = {cf.top_roughness * 1e9:.1f} nm',
                f'\nBottom roughness = {cf.bottom_roughness * 1e9:.1f} nm\n',
                f'\n{percentage:.0f}% of phonons reached the cold side\n'
        )
        file.writelines(info)


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

    retherm = 100*np.sum(scatter_stats.hot_side) / total
    internal = 100*np.sum(scatter_stats.internal) / total

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

    if cf.include_holes:
        scat_on_holes = 100*(np.sum(scatter_stats.hole_diffuse) +
                             np.sum(scatter_stats.hole_specular)) / total
        scat_on_holes_diff = 100*np.sum(scatter_stats.hole_diffuse) / total_hole
        scat_on_holes_spec = 100*np.sum(scatter_stats.hole_specular) / total_hole
        info2 = (
                f'\n{scat_on_holes:.2f}% - scattering on hole walls ',
                f'({scat_on_holes_diff:.2f}% - diffuse, ',
                f'{scat_on_holes_spec:.2f}% - specular)',
        )

    if cf.include_pillars:
        scat_on_pill = 100*(np.sum(scatter_stats.pillar_diffuse) +
                            np.sum(scatter_stats.pillar_specular)) / total
        scat_on_pill_diff = 100*np.sum(scatter_stats.pillar_diffuse) / total_pill
        scat_on_pill_spec = 100*np.sum(scatter_stats.pillar_specular) / total_pill
        info3 = (
                f'\n{scat_on_pill:.2f}% - scattering on pillar walls ',
                f'({scat_on_pill_diff:.2f}% - diffuse, ',
                f'{scat_on_pill_spec:.2f}% - specular)'
        )

    # Write info into a text file:
    with open("Information.txt", "a", encoding="utf-8") as file:
        file.writelines(info1)
        if cf.include_holes:
            file.writelines(info2)
        if cf.include_pillars:
            file.writelines(info3)
