"""Simple input file to simulate a 2D graphene sheet"""

OUTPUT_FOLDER_NAME             = "Graphene sheet"
NUMBER_OF_PHONONS              = 1000
NUMBER_OF_TIMESTEPS            = 60000
TIMESTEP                       = 2e-12
T                              = 300
THICKNESS                      = 150e-9
WIDTH                          = 1000e-9
LENGTH                         = 2200e-9

# Material:
MEDIA                          = "Graphite"
IS_TWO_DIMENSIONAL_MATERIAL    = True
SPECIFIC_HEAT_CAPACITY         = 710.0  # [J/kg/K] for Graphene at 300 K

PHONON_SOURCES = [Source(size_x=WIDTH, size_y=0, size_z=THICKNESS, angle_distribution="random_up")]
