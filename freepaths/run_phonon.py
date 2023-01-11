"""Module that runs one phonon through the structure"""


from parameters import *
from scattering import *
from events import ScatteringTypes


def run_phonon(phonon, flight, scatter_stats, segment_stats, thermal_maps, scatter_maps, material):
    """Run one phonon through the system and record parameters of this run"""

    scattering_types = ScatteringTypes()

    # Run the phonon step-by-step:
    for step_number in range(NUMBER_OF_TIMESTEPS):

        # If phonon has not reached the cold side yet:
        if phonon.is_in_system:

            # Check if different scattering events happened during current time step:
            if INCLUDE_INTERNAL_SCATTERING:
                internal_scattering(phonon, flight, scattering_types)
            surface_scattering(phonon, scattering_types)
            reinitialization(phonon, scattering_types)

            # Record scattering events if any:
            if scattering_types.is_scattered:
                scatter_stats.save_scattering_events(phonon.y, scattering_types)
                flight.add_point_to_path()

            # If no diffuse scattering event occurred, keep measuring the paths and time:
            if not (scattering_types.is_diffuse or scattering_types.is_internal):
                flight.add_step()

            # If diffuse scattering has occurred, reset phonon free path:
            else:
                flight.save_free_paths()
                flight.restart()
                phonon.assign_internal_scattering_time(material)

            # Update scattering and energy maps:
            if OUTPUT_SCATTERING_MAP and scattering_types.is_scattered:
                scatter_maps.add_scattering_to_map(phonon, scattering_types)
            thermal_maps.add_energy_to_maps(phonon, step_number, material)

            # Record time spent in the segment:
            segment_stats.record_time_in_segment(phonon.y)

            # Phonon makes a step forward:
            phonon.move()

            # Reset scattering types for the next step:
            scattering_types.reset()

        # If the phonon reached cold side, record a few parameters and break the loop:
        else:
            flight.add_point_to_path()
            flight.finish(step_number)
            break
