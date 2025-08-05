"""
Module that runs one particle through the structure, from particle source to a cold side.
After the run, we record various parameters of this run into flight object.
"""
from freepaths.config import cf
from freepaths.scattering import internal_scattering, surface_scattering, reinitialization
from freepaths.scattering_types import ScatteringTypes, ScatteringPlaces
from freepaths.particle_types import ParticleType

def run_particle(particle, flight, scatter_stats, places_stats, segment_stats, thermal_maps, scatter_maps, material):
    """Run one particle through the system and record parameters of this run"""

    # Initialize object that will store scattering types:
    scattering_types = ScatteringTypes()
    triangle_scattering_places = ScatteringPlaces()

    # Run the particle step-by-step:
    for step_number in range(cf.number_of_timesteps):

        # If the particle reached a cold side, record it and break the loop:
        if particle.has_crossed_cold_side:
            flight.add_point_to_path()
            flight.save_free_paths()
            flight.finish(step_number, cf.timestep)
            break

        # If the particle reached a hot side, record it and break the loop:
        if particle.has_crossed_hot_side and not cf.rethermalization_on_hot_sides:
            flight.add_point_to_path()
            flight.save_free_paths()
            break

        # Check if different scattering events happened during current time step:
        if cf.include_internal_scattering:
            internal_scattering(particle, flight, scattering_types)
        surface_scattering(particle, scattering_types, triangle_scattering_places)
        if cf.rethermalization_on_hot_sides:
            reinitialized = reinitialization(particle, scattering_types)
            if reinitialized:
                flight.reset_travel_time()

        # If any scattering has occurred, record it:
        if scattering_types.is_scattered:
            flight.add_point_to_path()
            scatter_stats.save_scattering_events(particle.y, scattering_types)
            if cf.output_scattering_map:
                scatter_maps.add_scattering_to_map(particle, scattering_types)

        # Otherwise, record only if animation is requested:
        else:
            if cf.output_path_animation:
                flight.add_point_to_path()

        # If diffuse scattering has occurred, reset particle free path:
        if scattering_types.is_diffuse or scattering_types.is_internal:
            flight.save_free_paths()
            flight.restart()
            particle.assign_internal_scattering_time(material)
            if cf.is_two_dimensional_material:
                particle.phi = 0.0

        # If hole scattering has occured, record it:
        if triangle_scattering_places.is_scattered:
            places_stats.save_scattering_events(particle.y, triangle_scattering_places)
        if scattering_types.is_diffuse_on_hole:
            flight.save_hole_diff_scattering_angle(particle.theta)
        if scattering_types.is_specular_on_hole:
            flight.save_hole_spec_scattering_angle(particle.theta)

        else:
            flight.add_step(cf.timestep)
        
        # Record presence of the particle at this timestep and move on:
        thermal_maps.add_energy_to_maps(particle, step_number, material)
        segment_stats.record_time_in_segment(particle.y)
        scattering_types.reset()
        triangle_scattering_places.reset()
        particle.move()

