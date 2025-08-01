"""This example show how FreePATHS can output animation of particle paths in the structure"""

# !! Animations are currently broken !!

# Basic parameters:
OUTPUT_FOLDER_NAME             = "Animated Paths"
NUMBER_OF_PARTICLES            = 50
NUMBER_OF_TIMESTEPS            = 300
TIMESTEP                       = 20e-12
T                              = 4
THICKNESS                      = 150e-9
WIDTH                          = 500e-9
LENGTH                         = 4000e-9

# Multiprocessing
NUMBER_OF_PROCESSES = 10

# Animation parameters:
OUTPUT_PATH_ANIMATION          = True
OUTPUT_TRAJECTORIES_OF_FIRST   = 30
OUTPUT_ANIMATION_FPS           = 24

# Particle source:
PARTICLE_SOURCES               = [Source(x=0, y=0, z=0, size_x=WIDTH,  size_y=0, size_z=THICKNESS, angle_distribution="random", angle=np.pi/2)]
