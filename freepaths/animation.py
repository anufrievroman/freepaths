"""Module that creates animations from recorded phonon paths"""

import sys
import os
import shutil
import imageio
import numpy as np
import matplotlib.pyplot as plt

from freepaths.config import cf


def generate_frames_xy():
    """Generate animation frames"""
    data = np.genfromtxt("Data/Phonon paths.csv",
                         unpack=False,
                         delimiter=',',
                         skip_header=1,
                         encoding='utf-8')

    # Create XY plots for each step of each phonon:
    frame_number = 0
    for phonon_num in range(cf.phonons_in_animation):
        x_coords = np.trim_zeros(data[:, 3 * phonon_num], trim='b')
        y_coords = np.trim_zeros(data[:, 3 * phonon_num + 1], trim='b')
        steps = np.shape(x_coords)[0]
        cmap = plt.get_cmap("tab10")
        for step in range(1, steps):
            fig, ax = plt.subplots()
            ax.plot(x_coords[:step], y_coords[:step], color=cmap(phonon_num), linewidth=0.5)
            ax.set_xlabel('X (μm)', fontsize=12)
            ax.set_ylabel('Y (μm)', fontsize=12)
            ax.set_xlim([-0.55*cf.width*1e6, 0.55*cf.width*1e6])
            ax.set_ylim([0, cf.length*1e6])
            ax.set_aspect('equal')
            fig.savefig(f"Frames/frame_{frame_number:0>6}.png", dpi=300, bbox_inches="tight")
            plt.close(fig)
            frame_number += 1


def generate_animation_xy():
    """Generate animation of phonon path in XY plane"""
    images = []
    filenames = os.listdir("Frames/")
    for filename in filenames:
        images.append(imageio.imread(f"Frames/{filename}"))
    imageio.mimsave("Animated paths XY.gif", images)


def delete_frames_xy():
    """Delete frames after creating the animation"""
    folder_path = "Frames/"
    shutil.rmtree(folder_path)
