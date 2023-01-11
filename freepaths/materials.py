"""Module that assigns physical properties according to chosen material"""

import numpy as np


class Material:
    """Material of the simulated media with certain physical properties"""

    def __init__(self, material):
        self.name = material

        N = 1000

        if self.name == 'Si':  # Ref. APL 95 161901 (2009)
            self.default_speed = 6000  # [m/s] This is the speed for Debye approximation
            self.density = 2330  # [kg/m^3]
            self.dispersion = np.zeros((N, 4))
            self.dispersion[:, 0] = [k * 12e9 / (N - 1) for k in range(N)]  # Wavevectors
            self.dispersion[:, 1] = [abs(1369.42 * k - 2.405e-8 * (k ** 2) - 9.70e-19 * (k ** 3)) for k in self.dispersion[:, 0]]  # LA branch
            self.dispersion[:, 2] = [abs(1081.74 * k - 7.711e-8 * (k ** 2) + 5.674e-19 * (k ** 3) + 7.967e-29 * (k ** 4)) for k in self.dispersion[:, 0]]  # TA branch
            self.dispersion[:, 3] = self.dispersion[:, 2]  # TA branch

        elif self.name == 'SiC':  # Ref. https://journals.aps.org/prb/pdf/10.1103/PhysRevB.50.17054
            self.default_speed = 6500  # [m/s] Need to change this probably!
            self.density = 3215  # [kg/m^3]
            self.dispersion = np.zeros((N, 4))
            self.dispersion[:, 0] = [k * 14414281503 / (N - 1) for k in range(N)]  # Wavevectors
            self.dispersion[:, 1] = [abs(-3.48834e-18 * (k ** 3) + 1.7604452e-08 * (k ** 2) + 1737.36296 * k) for k in self.dispersion[:, 0]]  # LA branch
            self.dispersion[:, 2] = [abs(-2.21696e-19 * (k ** 3) - 3.4366886e-08 * (k ** 2) + 1077.98941 * k) for k in self.dispersion[:, 0]]  # TA branch
            self.dispersion[:, 3] = self.dispersion[:, 2]  # TA branch

        elif self.name == 'Diamond':  # Ref. https://www.sciencedirect.com/science/article/pii/S0008622315003358
            self.default_speed = 20000  # [m/s]
            self.density = 3500  # [kg/m^3]
            self.dispersion = np.zeros((N, 4))
            self.dispersion[:, 0] = [k * 11707071561.7 / (N - 1) for k in range(N)]  # Wavevectors
            self.dispersion[:, 1] = [abs(-1.347265e-18 * (k ** 3) - 8.855338e-08 * (k ** 2) + 4309.95222 * k) for k in self.dispersion[:, 0]]  # LA branch
            self.dispersion[:, 2] = [abs(-5.042335e-18 * (k ** 3) - 4.104260e-08 * (k ** 2) + 3185.66561 * k) for k in self.dispersion[:, 0]]  # TA branch
            self.dispersion[:, 3] = self.dispersion[:, 2]  # TA branch

        elif self.name == 'AlN':  # Ref. Davydov, PRB 58, 12899 (1998)
            self.default_speed = 6200  # [m/s]
            self.density = 3255  # [kg/m^3]
            self.dispersion = np.zeros((N, 4))
            self.dispersion[:, 0] = [k * 12576399382.998995 / (N - 1) for k in range(N)]  # Wavevectors
            self.dispersion[:, 1] = [abs(-2.990977e-18 * (k ** 3) + 3.08258e-08 * (k ** 2) + 946.677 * k) for k in self.dispersion[:, 0]]  # LA branch
            self.dispersion[:, 2] = [abs(-3.047928e-18 * (k ** 3) - 2.08813e-08 * (k ** 2) + 1852.27 * k) for k in self.dispersion[:, 0]]  # TA branch
            self.dispersion[:, 3] = self.dispersion[:, 2]  # TA branch
        else:
            raise ValueError('Specified material does not exist in the database')
