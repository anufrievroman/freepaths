"""Pristine bulk Si — reproduces Priyadarshi et al. 2023 (J. Appl. Phys. 133, 054301) Figs. 3–4"""

from scipy.constants import electron_volt

# General parameters:
OUTPUT_FOLDER_NAME                  = 'Bulk_Si_pristine'
NUMBER_OF_PARTICLES                 = 10000
T                                   = 300.0
OUTPUT_SCATTERING_MAP               = False

# Simulation time parameters:
TIMESTEP                            = 1e-14
NUMBER_OF_TIMESTEPS                 = 2400000
NUMBER_OF_VIRTUAL_TIMESTEPS         = NUMBER_OF_TIMESTEPS * 4
NUMBER_OF_STABILIZATION_TIMEFRAMES  = 5
NUMBER_OF_TIMEFRAMES                = 8
LOW_MEMORY_USAGE                    = True

# Multiprocessing:
NUMBER_OF_PROCESSES                 = 10

# Electron parameters — match paper: 5 meV steps, 0–500 meV range, 15 nm MFP
ELECTRON_MFP                        = 15e-9       # [m] — pristine Si value from paper Fig. 7
MEAN_MAPPING_CONSTANT               = None        # auto-compute C from BTE calibration

# Material parameters:
MEDIA                               = 'Si'
MEDIA_FERMI_LEVEL                   = 0.1 * electron_volt  # 100 meV — paper's Fig. 3/4 example

# System dimensions [m]:
THICKNESS                           = 10000e-9
WIDTH                               = 10000e-9
LENGTH                              = 500e-9       # 500 nm channel, same as paper

# Map & profile parameters:
pixel_size = 50e-9
NUMBER_OF_PIXELS_X                  = int(WIDTH / pixel_size)
NUMBER_OF_PIXELS_Y                  = int(LENGTH / pixel_size)
IGNORE_FAULTY_PARTICLES             = False

# Particle source:
PARTICLE_SOURCES = [Source(x=0, y=0, z=0, size_x=WIDTH, size_y=0, size_z=THICKNESS, angle_distribution="random", angle=0)]
