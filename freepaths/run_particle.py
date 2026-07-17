"""
Module that runs one particle through the structure, from particle source to a cold side.
After the run, we record various parameters of this run into flight object.
"""
from freepaths.config import cf
from freepaths.scattering import internal_scattering, surface_scattering, reinitialization
from freepaths.scattering_types import ScatteringTypes, ScatteringPlaces
from freepaths.options import ParticleType, SimulationMode

def run_particle(particle, flight, scatter_stats, places_stats, segment_stats, thermal_maps, scatter_maps, material, mode):
    """Run one particle through the system and record parameters of this run"""

    # Initialize object that will store scattering types:
    scattering_types = ScatteringTypes()
    triangle_scattering_places = ScatteringPlaces()

    # Run the particle step-by-step:
    for step_number in range(cf.number_of_timesteps):

        # If the particle reached a cold side, record it and break the loop:
        if particle.has_crossed_cold_side:
            if not cf.low_memory_usage:
                flight.add_point_to_path()
            flight.save_free_paths()
            flight.finish(step_number, cf.timestep)
            break

        # If the particle reached a hot side, record it and break the loop:
        if particle.has_crossed_hot_side and not cf.rethermalization_on_hot_sides:
            if not cf.low_memory_usage:
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
            if not cf.low_memory_usage:
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
            # An inelastic (anharmonic) internal event rethermalizes the phonon: new
            # branch and frequency are drawn from the collision-rate-weighted
            # distribution. An elastic (impurity) internal event conserves the mode —
            # the direction randomization already applied above is all it does.
            # Applied only in the phonon tracing mode: in the MFP sampling mode each
            # phonon must keep its mode identity. Not applied in the gray approximation:
            if (scattering_types.is_internal
                    and cf.rethermalize_inelastic_scattering
                    and particle.type is ParticleType.PHONON
                    and mode is SimulationMode.PHONON_TRACING
                    and not cf.use_gray_approximation_mfp
                    and particle.internal_event_is_inelastic(material)):
                particle.rethermalize(material)
            particle.assign_internal_scattering_time(material)
            if cf.is_two_dimensional_material:
                particle.phi = 0.0

            # MFP sampling only: a short-tau, low-velocity phonon takes tiny hops and
            # essentially never reaches a boundary, but its mean free path already
            # converges after a modest number of scattering events, so stop early
            # rather than burning the full timestep budget on diminishing returns.
            # Tracing mode must keep going regardless, to build the flux/temperature map:
            if (mode is SimulationMode.PHONON_MFP_SAMPLING
                    and cf.max_number_of_scattering_events is not None
                    and len(flight.free_paths) >= cf.max_number_of_scattering_events):
                flight.finish(step_number, cf.timestep)
                break

        # If hole scattering has occurred, record it:
        if triangle_scattering_places.is_scattered:
            places_stats.save_scattering_events(particle.y, triangle_scattering_places)
        if scattering_types.is_diffuse_on_hole:
            flight.save_hole_diff_scattering_angle(particle.theta)
        if scattering_types.is_specular_on_hole:
            flight.save_hole_spec_scattering_angle(particle.theta)

        # If interface transmission has occurred, record it:
        if scattering_types.interfaces_transmission:
            places_stats.save_scattering_events(particle, scattering_types.interfaces_transmission)
        else:
            flight.add_step(cf.timestep)


        # Record presence of the particle at this timestep and move on:
        if thermal_maps is not None:
            thermal_maps.add_energy_to_maps(particle, step_number, material)
        segment_stats.record_time_in_segment(particle.y)
        scattering_types.reset()
        triangle_scattering_places.reset()
        particle.move()

