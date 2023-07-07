# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 17:27:20 2023

@author: Felix
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 17:52:36 2023

@author: Felix
"""

import numpy as np
from math import pi


# General parameters:
OUTPUT_FOLDER_NAME             = 'lattice left to right freepath'
NUMBER_OF_PHONONS              = 1000
NUMBER_OF_TIMESTEPS            = 60000
NUMBER_OF_NODES                = 400
TIMESTEP                       = 1/4*1e-12
T                              = 4.0
PLOTS_IN_TERMINAL              = True
OUTPUT_SCATTERING_MAP          = False
OUTPUT_RAW_THERMAL_MAP         = True
OUTPUT_TRAJECTORIES_OF_FIRST   = 100
NUMBER_OF_LENGTH_SEGMENTS      = 10


# Map & profiles parameters:
NUMBER_OF_PIXELS_X             = 200
NUMBER_OF_PIXELS_Y             = 200
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
WIDTH                          = 3120e-9
LENGTH                         = 7.8*400e-9




# Roughness [m]:
SIDE_WALL_ROUGHNESS            = 2e-9
HOLE_ROUGHNESS                 = 2e-9
TOP_ROUGHNESS                  = 0.2e-9
BOTTOM_ROUGHNESS               = 0.2e-9




# Hole array parameters [m]:
INCLUDE_HOLES                  = True
CIRCULAR_HOLE_DIAMETER         = 800e-9
INNER_CIRCULAR_HOLE_DIAMETER   = 800e-9 #1e-12 #300e-9
RECTANGULAR_HOLE_SIDE_X        = None#200e-9
RECTANGULAR_HOLE_SIDE_Y        = None#350e-9
ALPHA_ARC                      = 60*np.pi/180
ANGLE0                         = 60*np.pi/180
aa                             = ALPHA_ARC/6
d                              = 10e-9
PERIOD_X                       = (3.8+0.2)*CIRCULAR_HOLE_DIAMETER/2 
PERIOD_Y                       =  3.8*CIRCULAR_HOLE_DIAMETER/2 #+ d #3.3*np.sin(aa)*CIRCULAR_HOLE_DIAMETER
L_PERIOD_x                     = (3.8-1.8+0.2)*CIRCULAR_HOLE_DIAMETER/2

# Hot and cold sides [m]:
FREQUENCY_DETECTOR_SIZE        = WIDTH
FREQUENCY_DETECTOR_CENTER       = 0
sigma                          = 1e-9



INCLUDE_RIGHT_SIDEWALL         = False
INCLUDE_LEFT_SIDEWALL          = False
INCLUDE_TOP_SIDEWALL           = True
INCLUDE_BOTTOM_SIDEWALL        = True

# Hot and cold sides [m]:
COLD_SIDE_POSITION_RIGHT            = True
COLD_SIDE_POSITION_BOTTOM           = False
COLD_SIDE_POSITION_TOP              = False
HOT_SIDE_POSITION_LEFT              = True
HOT_SIDE_POSITION_BOTTOM            = False
PHONON_SOURCE_X                     = -WIDTH/2
PHONON_SOURCE_Y                     = 380e-9
PHONON_SOURCE_WIDTH_X               = 0
PHONON_SOURCE_WIDTH_Y               = 1e-9
PHONON_SOURCE_ANGLE_DISTRIBUTION    = 'random_right'

# Lattice of holes:
FIRST_HOLE_COORDINATE_X = -WIDTH/2-PERIOD_X #1.5/2*CIRCULAR_HOLE_DIAMETER
FIRST_HOLE_COORDINATE_Y =0# LENGTH/4
NUMBER_OF_PERIODS_X = 18
NUMBER_OF_PERIODS_Y = 3
HOLE_COORDINATES = np.zeros((NUMBER_OF_PERIODS_X * NUMBER_OF_PERIODS_Y, 3))
SCALING_FACTOR_RADIUS = np.zeros((NUMBER_OF_PERIODS_X * NUMBER_OF_PERIODS_Y, 1))
SCALING_FACTOR_INNER_RADIUS = np.zeros((NUMBER_OF_PERIODS_X * NUMBER_OF_PERIODS_Y, 1))
SCALE_ANGLE_V = 0.25
SCALE_ANGLE_H = 3.5/5
SCALE_ANGLE_H_REVERSE = 3.5/5
SCALE_ANGLE=[2,2.05,2.1,0.3,2.15,2.15]

HOLE_SHAPES = ['arccircle_v_lattice' for x in range(HOLE_COORDINATES.shape[0])]
hole_number = 0
for i in range(NUMBER_OF_PERIODS_Y):
    for j in range(NUMBER_OF_PERIODS_X):
        HOLE_COORDINATES[hole_number, 0] =  FIRST_HOLE_COORDINATE_X +j//6*PERIOD_X+ i*PERIOD_X/2#FIRST_HOLE_COORDINATE_X -(NUMBER_OF_PERIODS_X - 1) * PERIOD_X / 2 + j * PERIOD_X #-200*10**(-9)
        
        HOLE_COORDINATES[hole_number, 1] = FIRST_HOLE_COORDINATE_Y +i*PERIOD_Y  #+j//3*L_PERIOD_Y #FIRST_HOLE_COORDINATE + i * PERIOD_Y #1*10**(-6)
        HOLE_COORDINATES[hole_number, 2] = 0
        
        if hole_number% 6==0:
            SCALING_FACTOR_INNER_RADIUS[hole_number]= 2.2
            SCALING_FACTOR_RADIUS[hole_number]=2.6
        if hole_number%6==1:
            SCALING_FACTOR_INNER_RADIUS[hole_number]=2.8  #+sigma/LENGTH#/cos(aaa)
            SCALING_FACTOR_RADIUS[hole_number]=3.0
        if hole_number%6==2:
            SCALING_FACTOR_INNER_RADIUS[hole_number]=3.2  #+sigma/LENGTH#/cos(aaa)
            SCALING_FACTOR_RADIUS[hole_number]=3.4
        if hole_number%6 == 3:
            SCALING_FACTOR_INNER_RADIUS[hole_number]=3.6
            SCALING_FACTOR_RADIUS[hole_number]=3.8
        if hole_number% 6== 4:
             SCALING_FACTOR_INNER_RADIUS[hole_number]=3.6
             SCALING_FACTOR_RADIUS[hole_number]=3.8
             HOLE_SHAPES[hole_number] = 'arccircle_v_demi_down'
        if hole_number% 6== 5:
             SCALING_FACTOR_INNER_RADIUS[hole_number]=3.6
             SCALING_FACTOR_RADIUS[hole_number]=3.8
             HOLE_SHAPES[hole_number] = 'arccircle_v_demi_up'
            
        hole_number += 1

