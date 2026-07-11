from enum import Enum

class ParticleType(Enum):
    PHONON   = 0
    ELECTRON = 1

class SimulationMode(Enum):
    PHONON_TRACING      = 0
    PHONON_MFP_SAMPLING = 1
    ELECTRON            = 2