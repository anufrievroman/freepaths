# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 10:34:51 2023

@author: Felix
"""

import numpy as np
from math import pi


# General parameters:
OUTPUT_FOLDER_NAME             = 'multiple curve'
NUMBER_OF_PHONONS              = 300
NUMBER_OF_TIMESTEPS            = 30000
NUMBER_OF_NODES                = 400
TIMESTEP                       = 1/2*1e-12
T                              = 4.0
PLOTS_IN_TERMINAL              = True
OUTPUT_SCATTERING_MAP          = False
OUTPUT_RAW_THERMAL_MAP         = True
OUTPUT_TRAJECTORIES_OF_FIRST   = 300
NUMBER_OF_LENGTH_SEGMENTS      = 10



# Map & profiles parameters:
NUMBER_OF_PIXELS_X             = 400
NUMBER_OF_PIXELS_Y             = 400
NUMBER_OF_TIMEFRAMES           = 6


# Material parameters:
MEDIA                          = 'Si'
SPECIFIC_HEAT_CAPACITY         = 0.0176  # [J/kg/K] for Si at 4 K
#SPECIFIC_HEAT_CAPACITY        = 714  # [J/kg/K] for Si at 300 K
# SPECIFIC_HEAT_CAPACITY       = 606  # [J/kg/K] for SiC at 300 K


# Internal scattering:
INCLUDE_INTERNAL_SCATTERING    = False
USE_GRAY_APPROXIMATION_MFP     = False
GRAY_APPROXIMATION_MFP         = None


# System dimensions [m]:
THICKNESS                      = 150e-9
WIDTH                          = 800e-9
LENGTH                         = 800e-9




# Roughness [m]:
SIDE_WALL_ROUGHNESS            = 2e-9
HOLE_ROUGHNESS                 = 2e-9
TOP_ROUGHNESS                  = 0.2e-9
BOTTOM_ROUGHNESS               = 0.2e-9




# Hole array parameters [m]:
INCLUDE_HOLES                  = True
CIRCULAR_HOLE_DIAMETER         = 2*800e-9
INNER_CIRCULAR_HOLE_DIAMETER   = 2*800e-9 #1e-12 #300e-9
RECTANGULAR_HOLE_SIDE_X        = None#200e-9
RECTANGULAR_HOLE_SIDE_Y        = None#350e-9
ALPHA_ARC                      = 100*np.pi/180
ANGLE0                         = 25*np.pi/180

PERIOD_X                       = CIRCULAR_HOLE_DIAMETER/2
PERIOD_Y                       = CIRCULAR_HOLE_DIAMETER/2 #+ d #3.3*np.sin(aa)*CIRCULAR_HOLE_DIAMETER


# Hot and cold sides [m]:
e=94*1e-9
er=e/CIRCULAR_HOLE_DIAMETER*2
ar=1/3-er
FREQUENCY_DETECTOR_SIZE        = WIDTH
FREQUENCY_DETECTOR_CENTER       = 0
sigma                          = 1e-9



INCLUDE_RIGHT_SIDEWALL         = False
INCLUDE_LEFT_SIDEWALL          = False
INCLUDE_TOP_SIDEWALL           = False
INCLUDE_BOTTOM_SIDEWALL        = False

# Hot and cold sides [m]:
COLD_SIDE_POSITION_RIGHT            = True
COLD_SIDE_POSITION_LEFT             = True
COLD_SIDE_POSITION_TOP              = True
HOT_SIDE_POSITION_BOTTOM            = True
PHONON_SOURCE_X                     = -e/2
PHONON_SOURCE_Y                     = 0
PHONON_SOURCE_WIDTH_X               = e*0.9
PHONON_SOURCE_WIDTH_Y               = 0
PHONON_SOURCE_ANGLE_DISTRIBUTION    = 'random_up'

# Lattice of holes:
FIRST_HOLE_COORDINATE_X =-WIDTH/2-PERIOD_X/2 #1.5/2*CIRCULAR_HOLE_DIAMETER
FIRST_HOLE_COORDINATE_Y = 0
NUMBER_OF_PERIODS_X = 6
NUMBER_OF_PERIODS_Y = 2
HOLE_COORDINATES = np.zeros((NUMBER_OF_PERIODS_X * NUMBER_OF_PERIODS_Y, 3))
SCALING_FACTOR_RADIUS = np.zeros((NUMBER_OF_PERIODS_X * NUMBER_OF_PERIODS_Y, 1))
SCALING_FACTOR_INNER_RADIUS = np.zeros((NUMBER_OF_PERIODS_X * NUMBER_OF_PERIODS_Y, 1))

HOLE_SHAPES = ['arccircle_v_lattice_curve' for x in range(HOLE_COORDINATES.shape[0])]
hole_number = 0
for i in range(NUMBER_OF_PERIODS_Y):
    for j in range(NUMBER_OF_PERIODS_X):
        HOLE_COORDINATES[hole_number, 0] =  FIRST_HOLE_COORDINATE_X +j//3*PERIOD_X+ i*PERIOD_X/2#FIRST_HOLE_COORDINATE_X -(NUMBER_OF_PERIODS_X - 1) * PERIOD_X / 2 + j * PERIOD_X #-200*10**(-9)
        HOLE_COORDINATES[hole_number, 1] = FIRST_HOLE_COORDINATE_Y +i*PERIOD_Y  #+j//3*L_PERIOD_Y #FIRST_HOLE_COORDINATE + i * PERIOD_Y #1*10**(-6)
        HOLE_COORDINATES[hole_number, 2] = 0
        
        if hole_number% 3==0:
            SCALING_FACTOR_INNER_RADIUS[hole_number]= ar/2
            SCALING_FACTOR_RADIUS[hole_number]=ar
            HOLE_SHAPES[hole_number] = 'arccircle_v_lattice_curve_begin'
        if hole_number%3==1:
            SCALING_FACTOR_INNER_RADIUS[hole_number]=er+ar  
            SCALING_FACTOR_RADIUS[hole_number]=2*ar+er
        if hole_number%3==2:
            SCALING_FACTOR_INNER_RADIUS[hole_number]=2*(er+ar)  
            SCALING_FACTOR_RADIUS[hole_number]=3*ar+2*er
            
        hole_number += 1
