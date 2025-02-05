"""Module that creates animations from recorded phonon paths"""

import sys
import os
import shutil
import imageio
import numpy as np
import matplotlib.pyplot as plt

from freepaths.config import cf
from freepaths.output_structure import draw_structure_top_view


def generate_frames_xy():
    """Generate animation frames with phonon paths"""

    data = np.genfromtxt("Data/Phonon paths.csv",
                         unpack=False,
                         delimiter=',',
                         skip_header=1,
                         encoding='utf-8')

    # Create XY plots for each timestep where all phonon shown at the same time:
    number_of_steps = np.shape(data)[0]
    number_of_phonons = np.shape(data)[1]//3

    for step in range(1, number_of_steps):
        fig, ax = plt.subplots()

        # Draw the structure:
        patches = draw_structure_top_view(cf, color_back=cf.output_structure_color)

        for patch in patches:
            ax.add_patch(patch)

        # Draw the paths:
        for phonon_num in range(number_of_phonons):
            x_coords = np.trim_zeros(data[:, 3 * phonon_num], trim='b')
            y_coords = np.trim_zeros(data[:, 3 * phonon_num + 1], trim='b')
            ax.plot(x_coords[:step], y_coords[:step], linewidth=0.5)

        # Plot settings:
        ax.set_xlim([-0.55*cf.width*1e6, 0.55*cf.width*1e6])
        ax.set_ylim([0, cf.length*1e6])
        ax.set_aspect('equal')
        ax.axis('off')
        for spine in ['top', 'right', 'left', 'bottom']:
            ax.spines[spine].set_visible(False)
        fig.savefig(f"Frames/frame_{step:0>6}.png", dpi=600, bbox_inches="tight")
        plt.close(fig)

        # Progress:
        sys.stdout.write(f"\rAnimation: {step}/{number_of_steps - 1} frames")

    # Create XY plots for each step for each phonon one by one:
    # cmap = plt.get_cmap("tab10")
    # frame_number = 0
    # number_of_phonons = np.shape(data)[1]//3
    # for phonon_num in range(number_of_phonons):
        # x_coords = np.trim_zeros(data[:, 3 * phonon_num], trim='b')
        # y_coords = np.trim_zeros(data[:, 3 * phonon_num + 1], trim='b')
        # steps = np.shape(x_coords)[0]
        # for step in range(1, steps):
            # fig, ax = plt.subplots()
            # ax.plot(x_coords[:step], y_coords[:step], color=cmap(phonon_num), linewidth=0.5)
            # ax.axis('off')
            # for spine in ['top', 'right', 'left', 'bottom']:
                # ax.spines[spine].set_visible(False)
            # ax.set_xlim([-0.55*cf.width*1e6, 0.55*cf.width*1e6])
            # ax.set_ylim([0, cf.length*1e6])
            # ax.set_aspect('equal')
            # fig.savefig(f"Frames/frame_{frame_number:0>6}.png", dpi=600, bbox_inches="tight")
            # plt.close(fig)
            # frame_number += 1


def generate_animation_xy():
    """Generate animation of phonon path in XY plane"""
    sys.stdout.write(f"\rAnimation: creating animation file")
    images = []
    filenames = os.listdir("Frames/")
    for filename in filenames:
        images.append(imageio.imread(f"Frames/{filename}"))

    # Create a GIF file:
    imageio.mimsave("Animated paths XY.gif", images,
                    duration=1000 * 1/cf.output_animation_fps, subrectangles=True)

    # Create an MP4 file:
    # writer = imageio.get_writer('Animated paths XY.mp4', fps=cf.output_animation_fps)

    # for im in images:
        # writer.append_data(imageio.imread(im))
    # writer.close()


def delete_frames():
    """Delete frames after creating the animation"""
    folder_path = "Frames/"
    shutil.rmtree(folder_path)


def create_animation():
    """Main function that creates the animation"""
    generate_frames_xy()
    generate_animation_xy()
    # delete_frames()
