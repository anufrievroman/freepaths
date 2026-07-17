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

ELECTRON_MFP                     = 10e-9 # [m]
MEAN_MAPPING_CONSTANT            = 5e-6 # [m²]

FERMI_LEVEL_LOWER_BOUND          = -0.2  # [eV] lower end of post-processing Fermi level sweep
FERMI_LEVEL_UPPER_BOUND          =  0.1  # [eV] upper end of post-processing Fermi level sweep;
# far above the band edge the Fermi window leaves the sampled energy range and
# the MC results become noisy, so there is little point in sweeping further


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
USE_GRAY_APPROXIMATION_MFP       = False
GRAY_APPROXIMATION_MFP           = None

# Sample phonon frequencies and branches from the tabulated dispersion
# (weight k^2 dk per bin, heat capacity weighting, group velocity weighting at
# the source) instead of the legacy Debye-approximation rejection sampling.
# The particles are then equal-energy bundles and the thermal maps record a
# constant weight per phonon instead of h*w:
SAMPLE_FROM_DISPERSION           = True

# Re-draw phonon branch and frequency at inelastic (anharmonic) internal
# scattering events from the collision-rate-weighted distribution (Peraud &
# Hadjiconstantinou, PRB 84, 205331 (2011)). This restores local thermal
# equilibrium of the phonon population, which is required for correct
# Fourier-law thermal conductivity when internal scattering dominates (e.g.
# bulk-like structures at room temperature). Elastic internal events
# (impurity/alloy scattering) never rethermalize: they conserve the mode and
# only randomize the direction (see Material.phonon_scattering_rates).
# Applies only to the phonon tracing mode; MFP sampling and electron modes
# keep the particle identity by construction.
RETHERMALIZE_INELASTIC_SCATTERING = True

# Convert deposited particle energy into the temperature profile using the
# dispersion-only heat capacity (Material.dispersion_heat_capacity, summed only
# over the branches in the tabulated dispersion) instead of the experimental
# heat capacity fit (Material.assign_heat_capacity). On by default because it
# makes the reported temperature and kappa self-consistent with the
# dispersion-based sampling/rethermalization scheme and the model's own RTA
# integral. Set to False to use the real material heat capacity instead (e.g.
# including optical branches absent from the tabulated dispersion):
USE_DISPERSION_HEAT_CAPACITY     = True

# MFP sampling mode only (no effect in phonon tracing, where particles must keep
# going for the full budget to build the flux map): end a phonon's flight early
# once it has accumulated this many free-path segments, instead of always running
# to NUMBER_OF_TIMESTEPS. A short-tau, low-velocity phonon takes tiny hops and
# diffuses so slowly that it essentially never reaches a domain boundary, but its
# mean free path already converges after a modest number of scattering events, so
# running it to the full timestep budget is wasted compute. 1000 is on by default,
# since it was validated (bulk Si MFP-sampling convergence study, Data/BulkSi_Phonon_MFP/)
# to give the same result as 5000 while being far cheaper. Set to None to disable
# (run to the timestep budget or a boundary, as before):
MAX_NUMBER_OF_SCATTERING_EVENTS  = 1000

# Grain boundary scattering (polycrystalline materials):
# GRAIN_SIZE sets the mean grain diameter [m]. None disables grain boundary scattering entirely.
# GRAIN_SIZE_STD sets the standard deviation [m]; the grain size is drawn per phonon from a
# lognormal distribution with this mean and std. Set to 0 for monodisperse (single grain size).
# GRAIN_ROUGHNESS is the RMS disorder width of the grain boundary [m], used in a Soffer-type
# specularity factor: at low frequencies (long wavelengths) the boundary appears smooth and
# phonons pass through; at high frequencies it acts as a classical diffuse scatterer.
# Typical values: GRAIN_SIZE 100 nm – 10 µm, GRAIN_ROUGHNESS 0.5–2 nm.
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
