# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 16:28:41 2023

@author: Felix
"""



import numpy as np
from math import pi


# General parameters:
OUTPUT_FOLDER_NAME             = 'lattice down to up new'
NUMBER_OF_PHONONS              = 300
NUMBER_OF_TIMESTEPS            = 60000
NUMBER_OF_NODES                = 400
TIMESTEP                       = 1/4*1e-12
T                              = 4.0
PLOTS_IN_TERMINAL              = True
OUTPUT_SCATTERING_MAP          = False
OUTPUT_RAW_THERMAL_MAP         = True
OUTPUT_TRAJECTORIES_OF_FIRST   = 300
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
WIDTH                          = 800e-9
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
ALPHA_ARC                      = 100*np.pi/180
ANGLE0                         = 25*np.pi/180
aa                             = ALPHA_ARC/6
d                              = 10e-9
PERIOD_X                       = (1)*CIRCULAR_HOLE_DIAMETER/2 
PERIOD_Y                       = CIRCULAR_HOLE_DIAMETER/2 #+ d #3.3*np.sin(aa)*CIRCULAR_HOLE_DIAMETER
L_PERIOD_x                     = 0

# Hot and cold sides [m]:
FREQUENCY_DETECTOR_SIZE        = WIDTH
FREQUENCY_DETECTOR_CENTER       = 0
sigma                          = 1e-9



INCLUDE_RIGHT_SIDEWALL         = True
INCLUDE_LEFT_SIDEWALL          = True
INCLUDE_TOP_SIDEWALL           = False
INCLUDE_BOTTOM_SIDEWALL        = False

# Hot and cold sides [m]:
COLD_SIDE_POSITION_RIGHT            = False
COLD_SIDE_POSITION_LEFT             = False
COLD_SIDE_POSITION_TOP              = True
HOT_SIDE_POSITION_BOTTOM            = True
PHONON_SOURCE_X                     = -WIDTH/2+ 3.5*CIRCULAR_HOLE_DIAMETER/2
PHONON_SOURCE_Y                     = 0
PHONON_SOURCE_WIDTH_X               = 100e-9
PHONON_SOURCE_WIDTH_Y               = 0
PHONON_SOURCE_ANGLE_DISTRIBUTION    = 'random_up'

# Lattice of holes:
FIRST_HOLE_COORDINATE_X = -WIDTH/2#-PERIOD_X #1.5/2*CIRCULAR_HOLE_DIAMETER
FIRST_HOLE_COORDINATE_Y =  LENGTH/4
NUMBER_OF_PERIODS_X = 12
NUMBER_OF_PERIODS_Y = 2
HOLE_COORDINATES = np.zeros((NUMBER_OF_PERIODS_X * NUMBER_OF_PERIODS_Y, 3))
SCALING_FACTOR_RADIUS = np.zeros((NUMBER_OF_PERIODS_X * NUMBER_OF_PERIODS_Y, 1))
SCALING_FACTOR_INNER_RADIUS = np.zeros((NUMBER_OF_PERIODS_X * NUMBER_OF_PERIODS_Y, 1))
SCALE_ANGLE_V = 0.1
SCALE_ANGLE_H = 3.5/5
SCALE_ANGLE_H_REVERSE = 3.5/5
SCALE_ANGLE=[0.7,0.82,0.98,0.098,0.45,0.45]

HOLE_SHAPES = ['arccircle_v_lattice' for x in range(HOLE_COORDINATES.shape[0])]
hole_number = 0
for i in range(NUMBER_OF_PERIODS_Y):
    for j in range(NUMBER_OF_PERIODS_X):
        HOLE_COORDINATES[hole_number, 0] =  FIRST_HOLE_COORDINATE_X +j//6*PERIOD_X+ i*PERIOD_X/2#FIRST_HOLE_COORDINATE_X -(NUMBER_OF_PERIODS_X - 1) * PERIOD_X / 2 + j * PERIOD_X #-200*10**(-9)
        
        HOLE_COORDINATES[hole_number, 1] = FIRST_HOLE_COORDINATE_Y +i*PERIOD_Y  #+j//3*L_PERIOD_Y #FIRST_HOLE_COORDINATE + i * PERIOD_Y #1*10**(-6)
        HOLE_COORDINATES[hole_number, 2] = 0
        
        if hole_number% 3==0:
            SCALING_FACTOR_INNER_RADIUS[hole_number]= 1-2*94/800-2.90*172.7/800
            SCALING_FACTOR_RADIUS[hole_number]=1-94/800*2-2*172.7/800
        if hole_number%3==1:
            SCALING_FACTOR_INNER_RADIUS[hole_number]=1-94/800-2*172.7/800  #+sigma/LENGTH#/cos(aaa)
            SCALING_FACTOR_RADIUS[hole_number]=1-94/800-172.7/800
        if hole_number%3==2:
            SCALING_FACTOR_INNER_RADIUS[hole_number]=1-172.7/800  #+sigma/LENGTH#/cos(aaa)
            SCALING_FACTOR_RADIUS[hole_number]=1
        if hole_number%6 == 3:
            SCALING_FACTOR_INNER_RADIUS[hole_number]=1+94/800
            SCALING_FACTOR_RADIUS[hole_number]=1+94/800+127.7/800
        if hole_number% 6== 4:
             SCALING_FACTOR_INNER_RADIUS[hole_number]=1+94/800
             SCALING_FACTOR_RADIUS[hole_number]=1+94/800+0.9/10
             HOLE_SHAPES[hole_number] = 'arccircle_v_demi_up'
        if hole_number% 6== 5:
             SCALING_FACTOR_INNER_RADIUS[hole_number]=1+94/800
             SCALING_FACTOR_RADIUS[hole_number]=1+94/800+0.9/10
             HOLE_SHAPES[hole_number] = 'arccircle_v_demi_down'
            
        hole_number += 1

