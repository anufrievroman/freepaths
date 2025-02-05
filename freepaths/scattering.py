"""
Modules provides scattering processes on various complex objects.
Each function determines whether the scattering should happen and
call corresponding function for scattering on corresponding primitive.
"""

from math import sqrt

from freepaths.config import cf
from freepaths.move import move
from freepaths.scattering_primitives import *


def internal_scattering(ph, flight, scattering_types):
    """Check if the time passed since previous diffuse scattering event reached
    the time until an internal scattering event, and if yes, scatters randomly"""
    if flight.time_since_previous_scattering >= ph.time_of_internal_scattering:
        scattering_types.internal = random_scattering(ph)


def reinitialization(ph, scattering_types):
    """Re-thermalize (diffusely) phonon when it comes back to one of the hot sides"""
    x, y, _ = move(ph, cf.timestep)

    if cf.hot_side_position_bottom and y < 0:
        scattering_types.hot_side = horizontal_surface_up_scattering(ph, cf.side_wall_roughness, is_diffuse=True)

    if cf.hot_side_position_top and y > cf.length:
        scattering_types.hot_side = horizontal_surface_down_scattering(ph, cf.side_wall_roughness, is_diffuse=True)

    if cf.hot_side_position_right and x > cf.width / 2:
        scattering_types.hot_side = vertical_surface_left_scattering(ph, cf.side_wall_roughness, cf, is_diffuse=True)

    if cf.hot_side_position_left and x < -cf.width / 2:
        scattering_types.hot_side = vertical_surface_right_scattering(ph, cf.side_wall_roughness, cf, is_diffuse=True)


def scattering_on_right_sidewall(ph, scattering_types, x, y, z):
    """Scatter phonon if it reached right side wall"""
    if x > cf.width / 2:
        scattering_types.walls = vertical_surface_left_scattering(ph, cf.side_wall_roughness, cf)


def scattering_on_left_sidewall(ph, scattering_types, x, y, z):
    """Scatter phonon if it reached left side wall"""
    if x < -cf.width / 2:
        scattering_types.walls = vertical_surface_right_scattering(ph, cf.side_wall_roughness, cf)


def scattering_on_top_sidewall(ph, scattering_types, x, y, z):
    """Check if the phonon hits top side wall and output new vector"""
    if y > cf.length:
        scattering_types.walls = horizontal_surface_down_scattering(ph, cf.side_wall_roughness)


def scattering_on_bottom_sidewall(ph, scattering_types, x, y, z):
    """Check if the phonon hits bottom side wall and output new vector"""
    if y < 0.0:
        scattering_types.walls = horizontal_surface_up_scattering(ph, cf.side_wall_roughness)


def floor_scattering(ph, scattering_types, x, y, z):
    """Check if the phonon hits the floor surface and calculate new angles"""
    if z < -cf.thickness / 2:
        scattering_types.top_bottom = in_plane_surface_scattering(ph, cf.top_roughness)


def ceiling_scattering(ph, scattering_types, x, y, z):
    """Check if the phonon hits the ceiling surface and if this place has a pillar and output new vector"""
    if z <= cf.thickness / 2:
        return
    if cf.pillars:
        for pillar in cf.pillars:
            distance_from_pillar_center = sqrt(
                (x - pillar.x0) ** 2 + (y - pillar.y0) ** 2
            )
            is_under_pillar = distance_from_pillar_center < pillar.diameter / 2
            if is_under_pillar:
                if z > pillar.height + cf.thickness / 2:
                    scattering_types.top_bottom = in_plane_surface_scattering(ph, cf.top_roughness)
                return
        # Regular scattering if phonon is not under any of the pillars:
        if ph.z < cf.thickness / 2:
            scattering_types.top_bottom = in_plane_surface_scattering(ph, cf.top_roughness)

    else:
        scattering_types.top_bottom = in_plane_surface_scattering(ph, cf.top_roughness)


def surface_scattering(ph, scattering_types, triangle_scattering_places):
    """Check for a surface scattering on this step"""

    # Preliminary move to see if phonon would cross something:
    x, y, z = move(ph, cf.timestep)

    # Scattering on top and bottom surfaces:
    ceiling_scattering(ph, scattering_types, x, y, z)
    floor_scattering(ph, scattering_types, x, y, z)

    # Scattering on sidewalls:
    if cf.include_right_sidewall:
        scattering_on_right_sidewall(ph, scattering_types, x, y, z)
    if cf.include_left_sidewall:
        scattering_on_left_sidewall(ph, scattering_types, x, y, z)
    if cf.include_top_sidewall:
        scattering_on_top_sidewall(ph, scattering_types, x, y, z)
    if cf.include_bottom_sidewall:
        scattering_on_bottom_sidewall(ph, scattering_types, x, y, z)

    # Scattering on holes:
    if cf.holes:
        # Check for each hole and each hole type:
        for hole in cf.holes:
            if hole.is_inside(x, y, z, cf):
                hole.scatter(ph, scattering_types, x, y, z, cf)

            # If there was any scattering, then no need to check rest of the holes:
            if scattering_types.holes is not None:
                break

    # Check for each pillar:
    if cf.pillars:
        for pillar in cf.pillars:
            pillar.check_if_scattering(ph, scattering_types, x, y, z, cf)

            # If there was any scattering, then no need to check other pillars:
            if scattering_types.pillars is not None:
                break

    # Scattering on interfaces:
    if cf.interfaces:
        for interface in cf.interfaces:
            if interface.is_crossed(ph, x, y, z) and interface.is_transmitted():
                interface.scatter(ph, scattering_types, x, y, z, cf)

    # Correct angle if it became more than 180 degrees:
    ph.correct_angle()
