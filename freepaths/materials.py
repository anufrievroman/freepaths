"""
Module that assigns physical properties according to chosen material. References:
Si - Ref. APL 95 161901 (2009)
SiC - Ref. PRB 50 17054 (1994)
Diamond - Carbon 91 266-274 (2015)
AlN - Ref. PRB 58 12899 (1998)
Graphite/Graphene - Ref. Yanagisawa et al, Surface and Interface Analysis 37 133-136 (2005)
"""

import numpy as np
from freepaths.options import Materials


class Material:
    """Material of the simulated media with certain physical properties"""

    def __init__(self, material, num_points=1000):
        self.name = material

        if self.name == Materials.Si:
            A1 = 1369.42
            B1 = -2.405e-8
            C1 = -9.70e-19
            A2 = 1081.74
            B2 = -7.711e-8
            C2 = 5.674e-19
            D2 = 7.967e-29
            self.default_speed = 6000   # [m/s] This is the speed for Debye approximation
            self.density = 2330         # [kg/m^3]
            self.dispersion = np.zeros((num_points, 4))
            self.dispersion[:, 0] = [k * 12e9 / (num_points - 1) for k in range(num_points)]                         # Wavevectors
            self.dispersion[:, 1] = [abs(A1 * k + B1 * k**2 + C1 * k**3) for k in self.dispersion[:, 0]]             # LA branch
            self.dispersion[:, 2] = [abs(A2 * k + B2 * k**2 + C2 * k**3 + D2 * k**4) for k in self.dispersion[:, 0]] # TA branch
            self.dispersion[:, 3] = self.dispersion[:, 2]

        elif self.name == Materials.SiC:
            A1 = 1737.36296
            B1 = 1.7604452e-08
            C1 = -3.48834e-18
            A2 = 1077.98941
            B2 = -3.43668e-08
            C2 = -2.21696e-19
            self.default_speed = 6500   # [m/s] Need to change this probably!
            self.density = 3215         # [kg/m^3]
            self.dispersion = np.zeros((num_points, 4))
            self.dispersion[:, 0] = [k * 14414281503 / (num_points - 1) for k in range(num_points)]       # Wavevectors
            self.dispersion[:, 1] = [abs(C1 * k**3 + B1 * k**2 + A1 * k) for k in self.dispersion[:, 0]]  # LA branch
            self.dispersion[:, 2] = [abs(C2 * k**3 + B2 * k**2 + A2 * k) for k in self.dispersion[:, 0]]  # TA branch
            self.dispersion[:, 3] = self.dispersion[:, 2]

        elif self.name == Materials.Diamond:
            A1 = 4309.95222
            B1 = -8.855338e-08
            C1 = -1.347265e-18
            A2 = 3185.66561
            B2 = -4.104260e-08
            C2 = -5.042335e-18
            self.default_speed = 20000  # [m/s]
            self.density = 3500         # [kg/m^3]
            self.dispersion = np.zeros((num_points, 4))
            self.dispersion[:, 0] = [k * 11707071561.7 / (num_points - 1) for k in range(num_points)]     # Wavevectors
            self.dispersion[:, 1] = [abs(C1 * k**3 + B1 * k**2 + A1 * k) for k in self.dispersion[:, 0]]  # LA branch
            self.dispersion[:, 2] = [abs(C2 * k**3 + B2 * k**2 + A2 * k) for k in self.dispersion[:, 0]]  # TA branch
            self.dispersion[:, 3] = self.dispersion[:, 2]

        elif self.name == Materials.AlN:
            A1 = 946.677
            B1 = 3.08258e-08
            C1 = -2.990977e-18
            A2 = 1852.27
            B2 = -2.08813e-08
            C2 = -3.047928e-18
            self.default_speed = 6200     # [m/s]
            self.density = 3255           # [kg/m^3]
            self.dispersion = np.zeros((num_points, 4))
            self.dispersion[:, 0] = [k * 12576399382.998995 / (num_points - 1) for k in range(num_points)]  # Wavevectors
            self.dispersion[:, 1] = [abs(C1 * k**3 + B1 * k**2 + A1 * k) for k in self.dispersion[:, 0]]    # LA branch
            self.dispersion[:, 2] = [abs(C2 * k**3 + B2 * k**2 + A2 * k) for k in self.dispersion[:, 0]]    # TA branch
            self.dispersion[:, 3] = self.dispersion[:, 2]


        elif self.name == Materials.Graphite:
            A1 = 3640.918652525539
            B1 = -4.113040385347707e-08
            C1 = -1.2498941480986494e-18
            A2 = 2304.3679889997725
            B2 = -4.725351497510156e-08
            C2 = -1.522988691492939e-18
            A3 = 106.82509744986406
            B3 = 1.1143719974291039e-07
            C3 = -3.502550128712751e-18
            self.default_speed = 12900     # [m/s]
            self.density = 2230            # [kg/m^3]
            self.dispersion = np.zeros((num_points, 4))
            self.dispersion[:, 0] = [k * 14500000000 / (num_points - 1) for k in range(num_points)]         # Wavevectors
            self.dispersion[:, 1] = [abs(C1 * k**3 + B1 * k**2 + A1 * k) for k in self.dispersion[:, 0]]    # LA branch
            self.dispersion[:, 2] = [abs(C2 * k**3 + B2 * k**2 + A2 * k) for k in self.dispersion[:, 0]]    # TA branch
            self.dispersion[:, 3] = [abs(C3 * k**3 + B3 * k**2 + A3 * k) for k in self.dispersion[:, 0]]    # ZA branch

        else:
            raise ValueError('Specified material does not exist in the database')
