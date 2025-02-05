"""Default config file"""

import numpy as np
from freepaths.sources import Source


# General parameters:
OUTPUT_FOLDER_NAME               = "Si nanowire at 300 K"
NUMBER_OF_PHONONS                = 5000
NUMBER_OF_NODES                  = 400
T                                = 300
OUTPUT_SCATTERING_MAP            = False
OUTPUT_TRAJECTORIES_OF_FIRST     = 50
OUTPUT_STRUCTURE_COLOR           = "#F0F0F0"
NUMBER_OF_LENGTH_SEGMENTS        = 10

# Time parameters:
TIMESTEP                         = 1e-12
NUMBER_OF_TIMESTEPS              = 300000
NUMBER_OF_VIRTUAL_TIMESTEPS      = NUMBER_OF_TIMESTEPS * 3
NUMBER_OF_TIMEFRAMES             = 8
NUMBER_OF_STABILIZATION_TIMEFRAMES = 5


# Animation:
OUTPUT_PATH_ANIMATION            = False
OUTPUT_ANIMATION_FPS             = 24

# Map & profiles parameters:
NUMBER_OF_PIXELS_X               = 25
NUMBER_OF_PIXELS_Y               = 100
IGNORE_FAULTY_PHONONS            = False

# Material parameters:
MEDIA                            = "Si"

# Internal scattering:
INCLUDE_INTERNAL_SCATTERING      = True
USE_GRAY_APPROXIMATION_MFP       = False
GRAY_APPROXIMATION_MFP           = None

# System dimensions [m]:
THICKNESS                        = 150e-9
WIDTH                            = 500e-9
LENGTH                           = 2000e-9
IS_TWO_DIMENSIONAL_MATERIAL      = False
INCLUDE_RIGHT_SIDEWALL           = True
INCLUDE_LEFT_SIDEWALL            = True
INCLUDE_TOP_SIDEWALL             = False
INCLUDE_BOTTOM_SIDEWALL          = False

# Hot and cold sides [m]:
COLD_SIDE_POSITION_TOP           = True
COLD_SIDE_POSITION_BOTTOM        = False
COLD_SIDE_POSITION_RIGHT         = False
COLD_SIDE_POSITION_LEFT          = False
HOT_SIDE_POSITION_TOP            = False
HOT_SIDE_POSITION_BOTTOM         = True
HOT_SIDE_POSITION_RIGHT          = False
HOT_SIDE_POSITION_LEFT           = False

# Phonon source:
PHONON_SOURCES = [Source(x=0, y=0, z=0, size_x=0,  size_y=0, size_z=0, angle_distribution="random", angle=0)]

# Roughness [m]:
SIDE_WALL_ROUGHNESS              = 2e-9
HOLE_ROUGHNESS                   = 2e-9
PILLAR_ROUGHNESS                 = 2e-9
TOP_ROUGHNESS                    = 0.2e-9
BOTTOM_ROUGHNESS                 = 0.2e-9
PILLAR_TOP_ROUGHNESS             = 0.2e-9
INTERFACE_ROUGHNESS              = 0.2e-9

# Holes and pillars:
HOLES                            = []
PILLARS                          = []
INTERFACES                       = []

# Multiprocessing:
NUMBER_OF_PROCESSES              = 10
