"""Default config file"""

import numpy as np
from scipy.constants import k, electron_volt
from freepaths.sources import Source


# General parameters:
OUTPUT_FOLDER_NAME               = "Si nanowire at 300 K"
NUMBER_OF_PARTICLES              = 10000
NUMBER_OF_NODES                  = 400
T                                = 300
OUTPUT_SCATTERING_MAP            = False
OUTPUT_TRAJECTORIES_OF_FIRST     = 20
OUTPUT_STRUCTURE_COLOR           = "#F0F0F0"
NUMBER_OF_LENGTH_SEGMENTS        = 10
LOW_MEMORY_USAGE                 = False

# Time parameters:
TIMESTEP                         = 2e-12
NUMBER_OF_TIMESTEPS              = 200000
NUMBER_OF_VIRTUAL_TIMESTEPS      = NUMBER_OF_TIMESTEPS * 3
NUMBER_OF_TIMEFRAMES             = 8
NUMBER_OF_STABILIZATION_TIMEFRAMES = 5

# Electron parameters [eV]
IS_CARRIER_ELECTRON              = True
ENERGY_UPPER_BOUND               = 350e-3
ENERGY_LOWER_BOUND               = 0
ENERGY_STEP                      = 5e-3

ELECTRON_MFP                     = 10e-9 # [m] crystal (acoustic-phonon-limited) MFP, energy-independent
# Ionized-impurity scattering: if DOPING_CONCENTRATION > 0, the MFP becomes energy-
# and doping-dependent via Matthiessen's rule, 1/Lambda_tot(E) = 1/ELECTRON_MFP +
# 1/Lambda_ii(E), where Lambda_ii is the Brooks-Herring contribution from N_I ionized
# dopants (screened Coulomb scattering). This shortens the MFP with doping and raises
# the Seebeck coefficient. 0 -> pure crystal MFP.
DOPING_CONCENTRATION             = 0.0   # [m^-3] ionized-impurity (dopant) concentration N_I
# Surface depletion of carriers: with SURFACE_POTENTIAL > 0, a dead layer of width
# W = sqrt(2*eps_s*phi_s/(e*N_D)) (abrupt-junction depletion) is excluded around every
# free surface (walls, top/bottom, hole edges) for CARRIERS ONLY. phi_s is the surface
# band bending set by Fermi-level pinning at the native oxide / surface states. Needs
# a doping (DEPLETION_DOPING or DOPING_CONCENTRATION) and the material dielectric
# constant. 0 -> no depletion.
SURFACE_POTENTIAL                = 0.0   # [eV] surface band bending phi_s (0 = off)
# Doping used ONLY for the depletion-width formula. If 0/None it falls back to
# DOPING_CONCENTRATION. Set this (and leave DOPING_CONCENTRATION = 0) to add a surface
# dead layer from the real doping WITHOUT engaging Brooks-Herring, i.e. to isolate the
# depletion effect on top of a constant-MFP transport calibration.
DEPLETION_DOPING                 = 0.0   # [m^-3] dopant concentration for the dead-layer width only
# Legacy phenomenological energy dependence (only used when DOPING_CONCENTRATION = 0):
# Lambda(E) = ELECTRON_MFP * (E/kT)^p.  p=0 -> constant; p=2 -> ionized-impurity-like.
ELECTRON_MFP_ENERGY_EXPONENT     = 0.0
MEAN_MAPPING_CONSTANT            = 5e-6 # [m²]

FERMI_LEVEL_LOWER_BOUND          = -0.2  # [eV] lower end of post-processing Fermi level sweep
FERMI_LEVEL_UPPER_BOUND          =  0.1  # [eV] upper end of post-processing Fermi level sweep;

# Animation:
OUTPUT_PATH_ANIMATION            = False
OUTPUT_ANIMATION_FPS             = 24

# Map & profiles parameters:
NUMBER_OF_PIXELS_X               = 7
NUMBER_OF_PIXELS_Y               = 67
IGNORE_FAULTY_PARTICLES          = False
GRADIENT_FIT_RANGE               = (0.1, 0.9)

# Material parameters:
MEDIA                            = "Si"
MEDIA_FERMI_LEVEL                = None

# Internal scattering:
INCLUDE_INTERNAL_SCATTERING      = True

# Phonon hydrodynamic transport.
PHONON_HYDRODYNAMIC              = False
NUMBER_OF_HYDRODYNAMIC_PRERUNS   = 5
HYDRODYNAMIC_PRERUNS_WEIGHT      = 0.8
NUMBER_OF_HYDRODYNAMIC_PRERUN_PARTICLES = None

# Diagnostic control: when True, Normal events still fire but redraw the direction 
# randomly instead of biased toward the local drift
HYDRODYNAMIC_NORMAL_RESISTIVE    = False
ISOTOPE_C13_CONCENTRATION        = 0.0

MAX_NUMBER_OF_SCATTERING_EVENTS  = 1000

# Grain boundary scattering (polycrystalline materials):
GRAIN_SIZE                       = None
GRAIN_SIZE_STD                   = 0.0
GRAIN_ROUGHNESS                  = 1e-9

# System dimensions [m]:
THICKNESS                        = 150e-9
WIDTH                            = 200e-9
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
RETHERMALIZATION_ON_HOT_SIDES    = True

# Particle source:
PARTICLE_SOURCES = [Source(x=0, y=0, z=0, size_x=0,  size_y=0, size_z=0, angle_distribution="random", angle=0)]

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
BULKS                            = []

# Multiprocessing:
NUMBER_OF_PROCESSES              = 10
