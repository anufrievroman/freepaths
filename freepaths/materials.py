"""Module that assigns physical properties according to chosen material"""

import numpy as np


class Si:
    """Physical properties of silicon"""

    def __init__(self, num_points=1000):
        self.name = "Si"
        self.default_speed = 6000   # [m/s]
        self.density = 2330         # [kg/m^3]
        self.assign_dispersion(num_points)

    def assign_dispersion(self, num_points):
        """Assign phonon disperion. Ref. APL 95 161901 (2009)"""
        A1 = 1369.42
        B1 = -2.405e-8
        C1 = -9.70e-19
        A2 = 1081.74
        B2 = -7.711e-8
        C2 = 5.674e-19
        D2 = 7.967e-29
        self.dispersion = np.zeros((num_points, 4))
        self.dispersion[:, 0] = [k * 12e9 / (num_points - 1) for k in range(num_points)]      # Wavevectors
        self.dispersion[:, 1] = [abs(A1 * k + B1 * k**2 + C1 * k**3) for k in self.dispersion[:, 0]]             # LA branch
        self.dispersion[:, 2] = [abs(A2 * k + B2 * k**2 + C2 * k**3 + D2 * k**4) for k in self.dispersion[:, 0]] # TA branch
        self.dispersion[:, 3] = self.dispersion[:, 2]

    def relaxation_time(self, omega, temp):
        """Calculate relaxation time at a given frequency and temperature"""
        deb_temp = 152.0
        tau_impurity = 1 / (2.95e-45 * (omega ** 4))
        tau_umklapp = 1 / (0.95e-19 * (omega ** 2) * temp * np.exp(-deb_temp / temp))
        return 1 / ((1 / tau_impurity) + (1 / tau_umklapp))


class SiC:
    """Physical properties of silicon carbide"""

    def __init__(self, num_points=1000):
        self.name = "SiC"
        self.density = 3215         # [kg/m^3]
        self.default_speed = 6500   # [m/s] Need to change this probably!
        self.assign_dispersion(num_points)

    def assign_dispersion(self, num_points):
        """Assign phonon disperion. Ref. PRB 50 17054 (1994)"""
        A1 = 1737.36296
        B1 = 1.7604452e-08
        C1 = -3.48834e-18
        A2 = 1077.98941
        B2 = -3.43668e-08
        C2 = -2.21696e-19
        self.dispersion = np.zeros((num_points, 4))
        self.dispersion[:, 0] = [k * 14414281503 / (num_points - 1) for k in range(num_points)] # Wavevectors
        self.dispersion[:, 1] = [abs(C1 * k**3 + B1 * k**2 + A1 * k) for k in self.dispersion[:, 0]]  # LA branch
        self.dispersion[:, 2] = [abs(C2 * k**3 + B2 * k**2 + A2 * k) for k in self.dispersion[:, 0]]  # TA branch
        self.dispersion[:, 3] = self.dispersion[:, 2]


    def relaxation_time(self, omega, temp):
        """
        Calculate relaxation time at a given frequency and temperature.
        In SiC we also take into account 4 phonon scattering.
        Ref. Joshi et al, JAP 88, 265 (2000)
        """
        deb_temp = 1200
        tau_impurity = 1 / (8.46e-45 * (omega ** 4))
        tau_umklapp = 1 / (6.16e-20 * (omega ** 2) * temp * np.exp(-deb_temp / temp))
        tau_4p = 1 / (6.9e-23 * (temp ** 2) * (omega ** 2))
        return 1 / ((1 / tau_impurity) + (1 / tau_umklapp) + (1 / tau_4p))


class Graphite:
    """Physical properties of graphite"""

    def __init__(self, num_points=1000):
        self.name = "Graphite"
        self.density = 2230            # [kg/m^3]
        self.default_speed = 12900     # [m/s]
        self.assign_dispersion(num_points)

    def assign_dispersion(self, num_points):
        """Assign phonon disperion. Carbon 91 266-274 (2015)"""
        A1 = 3640.918652525539
        B1 = -4.113040385347707e-08
        C1 = -1.2498941480986494e-18
        A2 = 2304.3679889997725
        B2 = -4.725351497510156e-08
        C2 = -1.522988691492939e-18
        A3 = 106.82509744986406
        B3 = 1.1143719974291039e-07
        C3 = -3.502550128712751e-18
        self.dispersion = np.zeros((num_points, 4))
        self.dispersion[:, 0] = [k * 14500000000 / (num_points - 1) for k in range(num_points)] # Wavevectors
        self.dispersion[:, 1] = [abs(C1 * k**3 + B1 * k**2 + A1 * k) for k in self.dispersion[:, 0]] # LA branch
        self.dispersion[:, 2] = [abs(C2 * k**3 + B2 * k**2 + A2 * k) for k in self.dispersion[:, 0]] # TA branch
        self.dispersion[:, 3] = [abs(C3 * k**3 + B3 * k**2 + A3 * k) for k in self.dispersion[:, 0]] # ZA branch

    def relaxation_time(self, omega, temp):
        """
        Calculate relaxation time at a given frequency and temperature.
        In graphene, we assume that material is perfect, i.e. without impurity scattering.
        Ref. PRB 87, 115421 (2013)
        """
        deb_temp = 1000.0
        tau_umklapp = 1 / (3.18e-25 * (omega ** 2) * (temp ** 3) * np.exp(-deb_temp / (3*temp)))
        return 1 / ( 1 / tau_umklapp)


# Materials below are not fully supported and don't have the relaxation times:

class Diamond:
    """Physical properties of diamond"""

    def __init__(self, num_points=1000):
        self.name = "Diamond"
        self.density = 3500         # [kg/m^3]
        self.default_speed = 20000  # [m/s]
        self.assign_dispersion(num_points)

    def assign_dispersion(self, num_points):
        """Assign phonon disperion. Ref. PRB 58 12899 (1998)"""

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

    def relaxation_time(self, omega, temp):
        pass


class AlN:
    """Physical properties of AlN"""

    def __init__(self, num_points=1000):
        self.name = "AlN"
        self.density = 3255           # [kg/m^3]
        self.default_speed = 6200     # [m/s]
        self.dispersion = np.zeros((self.num_points, 4))
        self.assign_dispersion(num_points)

    def assign_dispersion(self, num_points):
        """Assign phonon disperion. Ref. Yanagisawa et al, Surface and Interface Analysis 37 133-136 (2005)"""
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

    def relaxation_time(self, omega, temp):
        pass

