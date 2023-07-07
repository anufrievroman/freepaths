# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 09:54:35 2023

@author: Felix
"""

# -*- coding: utf-8 -*-
"""
Created on Thu May 18 17:16:50 2023

@author: Felix
"""

import numpy as np
from math import pi , asin , cos , sin


# General parameters:
OUTPUT_FOLDER_NAME             = 'elbow_e_test_freepath'
NUMBER_OF_PHONONS              = 1000
NUMBER_OF_TIMESTEPS            = 30000
NUMBER_OF_NODES                = 400
TIMESTEP                       = 1e-12
T                              = 2.0
PLOTS_IN_TERMINAL              = True
OUTPUT_SCATTERING_MAP          = False
OUTPUT_RAW_THERMAL_MAP         = True
OUTPUT_TRAJECTORIES_OF_FIRST   = 100
NUMBER_OF_LENGTH_SEGMENTS      = 1000



# Map & profiles parameters:
NUMBER_OF_PIXELS_X             = 300
NUMBER_OF_PIXELS_Y             = 300
NUMBER_OF_TIMEFRAMES           = 6


# Material parameters:
MEDIA                          = 'Si'
SPECIFIC_HEAT_CAPACITY         = 0.0176  # [J/kg/K] for Si at 4 K
#SPECIFIC_HEAT_CAPACITY        = 714  # [J/kg/K] for Si at 300 K
# SPECIFIC_HEAT_CAPACITY       = 606  # [J/kg/K] for SiC at 300 K


# Internal scattering:
INCLUDE_INTERNAL_SCATTERING    = True
USE_GRAY_APPROXIMATION_MFP     = False
GRAY_APPROXIMATION_MFP         = None


# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 3000e-9
LENGTH                         = 600e-9 



# Roughness [m]:
SIDE_WALL_ROUGHNESS            = 2e-9
HOLE_ROUGHNESS                 = 2e-9
TOP_ROUGHNESS                  = 0.2e-9
BOTTOM_ROUGHNESS               = 0.2e-9




# Hole array parameters [m]:
INCLUDE_HOLES                  = True
CIRCULAR_HOLE_DIAMETER         = 4e-9
RECTANGULAR_HOLE_SIDE_X        = 3000e-9
RECTANGULAR_HOLE_SIDE_Y        = 500e-9
e                              = 220e-9
h                              = e
PERIOD_X                       = h + RECTANGULAR_HOLE_SIDE_X
PERIOD_Y                       = 0




# Hot and cold sides [m]:
sigma                          = 200e-9
FREQUENCY_DETECTOR_SIZE        = h #3*WIDTH
FREQUENCY_DETECTOR_CENTER      = 0



INCLUDE_RIGHT_SIDEWALL         = False
INCLUDE_LEFT_SIDEWALL          = False
INCLUDE_TOP_SIDEWALL           = False
INCLUDE_BOTTOM_SIDEWALL        = True

# Hot and cold sides [m]:
COLD_SIDE_POSITION_RIGHT            = True
COLD_SIDE_POSITION_TOP              = True
HOT_SIDE_POSITION_LEFT              = True
HOT_SIDE_POSITION_BOTTOM            = False
PHONON_SOURCE_X                     = -WIDTH/2
PHONON_SOURCE_Y                     =  e/2
PHONON_SOURCE_WIDTH_X               = 0
PHONON_SOURCE_WIDTH_Y               = e
PHONON_SOURCE_ANGLE_DISTRIBUTION    = 'random_right'

# Lattice of holes:
FIRST_HOLE_COORDINATE_X = -RECTANGULAR_HOLE_SIDE_X /2-h/2
FIRST_HOLE_COORDINATE_Y = e + RECTANGULAR_HOLE_SIDE_Y/2
NUMBER_OF_PERIODS_X = 2
NUMBER_OF_PERIODS_Y = 1
HOLE_COORDINATES = np.zeros((2, 3))

HOLE_SHAPES = ['rectangle','rectangle']


hole_number = 0
for i in range(NUMBER_OF_PERIODS_Y):
    for j in range(NUMBER_OF_PERIODS_X):
        HOLE_COORDINATES[hole_number, 0] =  FIRST_HOLE_COORDINATE_X +j*PERIOD_X
        HOLE_COORDINATES[hole_number, 1] = FIRST_HOLE_COORDINATE_Y   #+j//3*L_PERIOD_Y #FIRST_HOLE_COORDINATE + i * PERIOD_Y #1*10**(-6)
        HOLE_COORDINATES[hole_number, 2] = 0
        
        
            
        hole_number += 1