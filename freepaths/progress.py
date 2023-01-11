"""Module that displays progress bar"""

import sys
from parameters import *


class Progress:
    """Progress bar for the program run"""

    def __init__(self):
        self.percentage = -1

    def render(self, phonon_number):
        """Check if percentage increased since last iteration and display the progress"""
        new_percentage = 100 * phonon_number // NUMBER_OF_PHONONS
        if new_percentage > self.percentage:
            self.percentage = new_percentage
            sys.stdout.write('\r' + 'Progress: ' + str(self.percentage) + '%')
            sys.stdout.flush()
