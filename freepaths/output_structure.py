"""Module that draws shape of the structure"""

from matplotlib.patches import Rectangle, Circle, Polygon
from freepaths.scatterers import *


def draw_structure(cf, color_holes="white", color_back="gray"):
    """Draw shape of the structure using patches from matplotlib"""

    # Overal shape as gray area:
    patches = [Rectangle((-1e6*cf.width/2, 0), 1e6*cf.width, 1e6*cf.length, facecolor=color_back)]

    # Holes as white patches:
    for hole in cf.holes:

        # Circular hole:
        if isinstance(hole, CircularHole):
            patch = Circle((1e6 * hole.x, 1e6 * hole.y), 1e6 * hole.diameter/2, facecolor=color_holes)
            patches.append(patch)

        # Rectangular hole:
        elif isinstance(hole, RectangularHole):
            patch = Rectangle((1e6 * (hole.x - hole.size_x/2), 1e6 * (hole.y - hole.size_y/2)),
                              1e6 * hole.size_x, 1e6 * hole.size_y, facecolor=color_holes)
            patches.append(patch)

        # Triangular hole up:
        elif isinstance(hole, TriangularUpHole):
            patch = Polygon([[1e6 * (hole.x - hole.size_x/2), 1e6 * (hole.y - hole.size_y/2)],
                             [1e6 * (hole.x + hole.size_x/2), 1e6 * (hole.y - hole.size_y/2)],
                             [1e6 * hole.x, 1e6 * (hole.y + hole.size_y/2)]], closed=True, facecolor=color_holes)
            patches.append(patch)

        # Triangular hole down:
        elif isinstance(hole, TriangularDownHole):
            patch = Polygon([[1e6 * (hole.x - hole.size_x/2), 1e6 * (hole.y + hole.size_y/2)],
                             [1e6 * (hole.x + hole.size_x/2), 1e6 * (hole.y + hole.size_y/2)],
                             [1e6 * hole.x, 1e6 * (hole.y - hole.size_y/2)]], closed=True, facecolor=color_holes)
            patches.append(patch)

        # Triangular half hole up (right):
        elif isinstance(hole, TriangularUpHalfHole) and hole.right_half:
            patch = Polygon([[1e6 * (hole.x), 1e6 * (hole.y - hole.size_y/2)],
                             [1e6 * (hole.x + hole.size_x/2), 1e6 * (hole.y - hole.size_y/2)],
                             [1e6 * hole.x, 1e6 * (hole.y + hole.size_y/2)]], closed=True, facecolor=color_holes)
            patches.append(patch)

        # Triangular half hole up (left):
        elif isinstance(hole, TriangularUpHalfHole) and not hole.right_half:
            patch = Polygon([[1e6 * (hole.x - hole.size_x/2), 1e6 * (hole.y - hole.size_y/2)],
                             [1e6 * (hole.x), 1e6 * (hole.y - hole.size_y/2)],
                             [1e6 * hole.x, 1e6 * (hole.y + hole.size_y/2)]], closed=True, facecolor=color_holes)
            patches.append(patch)

        # Triangular half hole down (right):
        elif isinstance(hole, TriangularDownHalfHole) and hole.right_half:
            patch = Polygon([[1e6 * (hole.x), 1e6 * (hole.y + hole.size_y/2)],
                             [1e6 * (hole.x + hole.size_x/2), 1e6 * (hole.y + hole.size_y/2)],
                             [1e6 * hole.x, 1e6 * (hole.y - hole.size_y/2)]], closed=True, facecolor=color_holes)
            patches.append(patch)

        # Triangular half hole down (left):
        elif isinstance(hole, TriangularDownHalfHole) and not hole.right_half:
            patch = Polygon([[1e6 * (hole.x - hole.size_x/2), 1e6 * (hole.y + hole.size_y/2)],
                             [1e6 * (hole.x), 1e6 * (hole.y + hole.size_y/2)],
                             [1e6 * hole.x, 1e6 * (hole.y - hole.size_y/2)]], closed=True, facecolor=color_holes)
            patches.append(patch)

    # Pillars as white patches:
    for pillar in cf.pillars:

        # Circular pillar:
        if isinstance(pillar, CircularPillar):
            patch = Circle((1e6 * pillar.x, 1e6 * pillar.y), 1e6 * pillar.diameter/2, facecolor=color_holes)
            patches.append(patch)

    # Phonon source areas as red patches:
    for source in cf.phonon_sources:
        width_x = 1e6*cf.width/49 if source.size_x == 0 else 1e6*source.size_x
        width_y = 1e6*cf.width/50 if source.size_y == 0 else 1e6*source.size_y
        x = 1e6*source.x - width_x / 2
        y = 1e6*source.y - width_y / 2
        phonon_source_patch = Rectangle((x, y), width_x, width_y, facecolor='red', alpha=0.5)
        patches.append(phonon_source_patch)

    return patches
