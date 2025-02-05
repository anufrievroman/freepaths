"""Module that draws shape of the structure"""

from matplotlib.patches import Rectangle
from freepaths.config import cf
from freepaths.scatterers import HorizontalPlane, VerticalPlane

def draw_structure_top_view(cf, color_holes="white", color_back="gray"):
    """Draw shape of the structure using patches from matplotlib"""

    # Overal shape as gray area:
    patches = [
        Rectangle(
            (-1e6 * cf.width / 2, 0),
            1e6 * cf.width,
            1e6 * cf.length,
            facecolor=color_back,
        )
    ]

    # Holes as white patches:
    for hole in cf.holes:
        patch = hole.get_patch(color_holes, cf)
        patches.extend(patch if isinstance(patch, list) else [patch])

    # Pillars as white patches:
    for pillar in cf.pillars:
        patch = pillar.get_patch(color_holes, cf)
        patches.extend(patch if isinstance(patch, list) else [patch])

    # Interfaces as white patches:
    for interface in cf.interfaces:
        if isinstance(interface, VerticalPlane):
            patch = interface.get_patch(color_holes, cf)
            patches.extend(patch if isinstance(patch, list) else [patch])

    # Phonon source areas as red patches:
    for source in cf.phonon_sources:
        width_x = 1e6 * cf.width / 49 if source.size_x == 0 else 1e6 * source.size_x
        width_y = 1e6 * cf.width / 50 if source.size_y == 0 else 1e6 * source.size_y
        x = 1e6 * source.x - width_x / 2
        y = 1e6 * source.y - width_y / 2
        phonon_source_patch = Rectangle(
            (x, y), width_x, width_y, facecolor='red', alpha=0.5
        )
        patches.append(phonon_source_patch)

    return patches


def draw_structure_side_view(cf, color_holes="white", color_back="gray"):
    """Draw shape of the structure using patches from matplotlib"""

    # Overal shape as gray area:
    patches = [
        Rectangle(
            (0, -1e6 * cf.thickness/2),
            1e6 * cf.length,
            1e6 * cf.thickness,
            facecolor=color_back,
        )
    ]

    # Interfaces as white patches:
    for interface in cf.interfaces:
        if isinstance(interface, HorizontalPlane):
            patch = interface.get_patch(color_holes, cf)
            patches.extend(patch if isinstance(patch, list) else [patch])

    return patches
