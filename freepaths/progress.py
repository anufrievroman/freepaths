"""Module that displays progress bar"""

import sys


class Progress:
    """Progress bar for the program run"""

    def __init__(self):
        self.percentage = -1

    def render(self, particle_number, number_of_particles):
        """Check if percentage increased since last iteration and display the progress"""
        new_percentage = 100 * particle_number // number_of_particles
        if new_percentage > self.percentage:
            self.percentage = new_percentage
            sys.stdout.write('\r' + 'Progress: ' + str(self.percentage) + '%')
            sys.stdout.flush()
