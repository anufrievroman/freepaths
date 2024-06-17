"""Module that assigns physical properties according to chosen material"""

import numpy as np


class Si:
    """
    Physical properties of silicon.
    Dispersion - Ref. APL 95 161901 (2009)
    Relaxation time - Maire et al, Scientific Reports 7, 41794 (2017)
    Heat capacity - Desai P.D. Journal of Physical and Chemical Reference Data 15, 67 (1986)
    """

    def __init__(self, temp, num_points=1000):
        self.name = "Si"
        self.default_speed = 6000   # [m/s]
        self.density = 2330         # [kg/m^3]
        self.temp = temp
        self.assign_dispersion(num_points)
        self.assign_heat_capacity()

    def assign_dispersion(self, num_points):
        """Assign phonon dispersion"""

        coefficients_LA = [-9.70e-19, -2.405e-8, 1369.42, 0]
        coefficients_TA = [7.967e-29, 5.674e-19, -7.711e-8, 1081.74, 0]

        self.dispersion = np.zeros((num_points, 4))
        self.dispersion[:, 0] = np.linspace(0, 12e9, num_points)  # Wavevectors
        self.dispersion[:, 1] = np.abs(np.polyval(coefficients_LA, self.dispersion[:, 0]))  # LA branch
        self.dispersion[:, 2] = np.abs(np.polyval(coefficients_TA, self.dispersion[:, 0]))  # TA branch
        self.dispersion[:, 3] = self.dispersion[:, 2]

    def relaxation_time(self, omega):
        """Calculate relaxation time at a given frequency and temperature"""
        deb_temp = 152.0
        tau_impurity = 1 / (2.95e-45 * (omega ** 4))
        tau_umklapp = 1 / (0.95e-19 * (omega ** 2) * self.temp * np.exp(-deb_temp / self.temp))
        return 1 / ((1 / tau_impurity) + (1 / tau_umklapp))

    def assign_heat_capacity(self):
        """Calculate heat capacity [J/kg/K] in 3 - 300K range using the polynomial fits"""
        below_20K_coeffs = np.array([0.00044801, -0.00239681,  0.00756769])
        between_20_and_50K_coeffs = np.array([-9.26222400e-04, 1.49879304e-01, -4.37458293e+00, 3.84245589e+01])
        above_50K_coeffs = np.array([-2.75839317e-06, -5.16662077e-03, 4.66701391e+00, -1.49876958e+02])
        if self.temp < 20:
            coeffs = below_20K_coeffs
        elif 20 <= self.temp <= 50:
            coeffs = between_20_and_50K_coeffs
        else:
            coeffs = above_50K_coeffs
        self.heat_capacity = np.polyval(coeffs, self.temp)


class SiC:
    """
    Physical properties of silicon carbide
    Dispersion - PRB 50 17054 (1994)
    Relaxation time - Joshi et al, JAP 88, 265 (2000)
    Heat capacity - Collins et al. Journal of Applied Physics 68, 6510 (1990)
    """

    def __init__(self, temp, num_points=1000):
        self.name = "SiC"
        self.density = 3215         # [kg/m^3]
        self.default_speed = 6500   # [m/s] Need to change this probably...
        self.temp = temp
        self.assign_dispersion(num_points)
        self.assign_heat_capacity()

    def assign_dispersion(self, num_points):
        """Assign phonon dispersion"""

        coefficients_LA = [-3.48834e-18, 1.7604452e-08, 1737.36296, 0]
        coefficients_TA = [-2.21696e-19, -3.43668e-08, 1077.98941, 0]

        self.dispersion = np.zeros((num_points, 4))
        self.dispersion[:, 0] = np.linspace(0, 14414281503, num_points)  # Wavevectors
        self.dispersion[:, 1] = np.abs(np.polyval(coefficients_LA, self.dispersion[:, 0]))  # LA branch
        self.dispersion[:, 2] = np.abs(np.polyval(coefficients_TA, self.dispersion[:, 0]))  # TA branch
        self.dispersion[:, 3] = self.dispersion[:, 2]

    def relaxation_time(self, omega):
        """Calculate relaxation time at a given frequency and temperature including 4 phonon scattering"""
        deb_temp = 1200
        tau_impurity = 1 / (8.46e-45 * (omega ** 4))
        tau_umklapp = 1 / (6.16e-20 * (omega ** 2) * self.temp * np.exp(-deb_temp / self.temp))
        tau_4p = 1 / (6.9e-23 * (self.temp ** 2) * (omega ** 2))
        return 1 / ((1 / tau_impurity) + (1 / tau_umklapp) + (1 / tau_4p))


    def assign_heat_capacity(self):
        """Calculate heat capacity [J/kg/K] in 3 - 500K range using the polynomial fits of experimental data"""

        below_90K_coeffs = np.array([7.08396921e-05, -9.66654246e-04, 9.03926727e-02, 2.99362037e-01])
        between_90_and_200K_coeffs = np.array([-6.59772636e-04, 3.02766713e-01, -4.12089642e+01, 1.82434354e+03])
        above_200K_coeffs = np.array([-3.57059689e-06, 1.21917876e-03, 2.28676930e+00, -7.23941447e+01])
        if self.temp < 90:
            coeffs = below_90K_coeffs
        elif 90 <= self.temp <= 200:
            coeffs = between_90_and_200K_coeffs
        else:
            coeffs = above_200K_coeffs
        self.heat_capacity = np.polyval(coeffs, self.temp)


class Graphite:
    """
    Physical properties of graphite.
    Dispersion - Carbon 91 266-274 (2015)
    Relaxation time - Ref. PRB 87, 115421 (2013)
    Heat capacity - Isaacs, L.L.; Wang, W.Y., Therm. Conduct. 17th, 55-61 (1981)
    """

    def __init__(self, temp, num_points=1000):
        self.name = "Graphite"
        self.density = 2230            # [kg/m^3]
        self.default_speed = 12900     # [m/s]
        self.temp = temp
        self.assign_dispersion(num_points)
        self.assign_heat_capacity()

    def assign_dispersion(self, num_points):
        """Assign phonon dispersion"""

        coefficients_LA = [-1.24989e-18, -4.11304e-08, 3640.918, 0]
        coefficients_TA = [-1.52298e-18, -4.72535e-08, 2304.367, 0]
        coefficients_ZA = [-3.50255e-18, 1.114371e-07, 106.8250, 0]

        self.dispersion = np.zeros((num_points, 4))
        self.dispersion[:, 0] = np.linspace(0, 14500000000, num_points)  # Wavevectors
        self.dispersion[:, 1] = np.abs(np.polyval(coefficients_LA, self.dispersion[:, 0]))  # LA branch
        self.dispersion[:, 2] = np.abs(np.polyval(coefficients_TA, self.dispersion[:, 0]))  # TA branch
        self.dispersion[:, 3] = np.abs(np.polyval(coefficients_ZA, self.dispersion[:, 0]))  # ZA branch

    def relaxation_time(self, omega):
        """
        Calculate relaxation time at a given frequency and temperature.
        In graphite, we assume that material is perfect, i.e. without impurity scattering.
        """
        deb_temp = 1000.0
        tau_umklapp = 1 / (3.18e-25 * (omega ** 2) * (self.temp ** 3) * np.exp(-deb_temp / (3*self.temp)))
        return 1 / ( 1 / tau_umklapp)


    def assign_heat_capacity(self):
        """Calculate heat capacity [J/kg/K] from the equation"""
        coeffs = np.array([6.309e-9, 6.27e-6, 8.729e-4, 0])
        self.heat_capacity = 1000 * np.polyval(coeffs, self.temp)


# Materials below are not fully supported and don't have the relaxation times:

class Diamond:
    """
    Physical properties of diamond
    Dispersion - Ref. PRB 58 12899 (1998)
    """

    def __init__(self, temp, num_points=1000):
        self.name = "Diamond"
        self.density = 3500         # [kg/m^3]
        self.default_speed = 20000  # [m/s]
        self.temp = temp
        self.assign_dispersion(num_points)

    def assign_dispersion(self, num_points):
        """Assign phonon dispersion"""

        A1 = 4309.95222
        B1 = -8.855338e-08
        C1 = -1.347265e-18
        A2 = 3185.66561
        B2 = -4.104260e-08
        C2 = -5.042335e-18

        self.dispersion = np.zeros((num_points, 4))
        self.dispersion[:, 0] = [k * 11707071561.7 / (num_points - 1) for k in range(num_points)]     # Wavevectors
        self.dispersion[:, 1] = [abs(C1 * k**3 + B1 * k**2 + A1 * k) for k in self.dispersion[:, 0]]  # LA branch
        self.dispersion[:, 2] = [abs(C2 * k**3 + B2 * k**2 + A2 * k) for k in self.dispersion[:, 0]]  # TA branch
        self.dispersion[:, 3] = self.dispersion[:, 2]

    def relaxation_time(self, omega):
        pass


class AlN:
    """
    Physical properties of AlN.
    Dispersion - Yanagisawa et al, Surface and Interface Analysis 37 133-136 (2005)
    """

    def __init__(self, temp, num_points=1000):
        self.name = "AlN"
        self.density = 3255           # [kg/m^3]
        self.default_speed = 6200     # [m/s]
        self.temp = temp              # [K]
        self.dispersion = np.zeros((self.num_points, 4))
        self.assign_dispersion(num_points)

    def assign_dispersion(self, num_points):
        """Assign phonon dispersion"""

        A1 = 946.677
        B1 = 3.08258e-08
        C1 = -2.990977e-18
        A2 = 1852.27
        B2 = -2.08813e-08
        C2 = -3.047928e-18

        self.dispersion[:, 0] = [k * 12576399382.998995 / (num_points - 1) for k in range(num_points)]
        self.dispersion[:, 1] = [abs(C1 * k**3 + B1 * k**2 + A1 * k) for k in self.dispersion[:, 0]]    # LA branch
        self.dispersion[:, 2] = [abs(C2 * k**3 + B2 * k**2 + A2 * k) for k in self.dispersion[:, 0]]    # TA branch
        self.dispersion[:, 3] = self.dispersion[:, 2]

    def relaxation_time(self, omega):
        pass

