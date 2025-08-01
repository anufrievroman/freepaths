from enum import Enum

class ParticleType(str, Enum):
    ELECTRON = "electron"
    PHONON   = "phonon"