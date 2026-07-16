"""Module that assigns physical properties according to chosen material"""

from abc import ABC, abstractmethod
import numpy as np
from scipy.constants import electron_volt, electron_mass, k as k_B, hbar, pi


class Material(ABC):

    # Names of the branches in the dispersion table, in column order:
    dispersion_branch_names = ['LA', 'TA1', 'TA2']

    @abstractmethod
    def assign_phonon_dispersion(self, num_points):
        """Assign phonon dispersion"""
        pass

    @abstractmethod
    def phonon_relaxation_time(self, omega):
        """Calculate relaxation time at a given frequency and temperature"""
        pass

    def phonon_scattering_rates(self, omega):
        """
        Return (inelastic_rate, elastic_rate) [1/s] at angular frequency omega.
        Inelastic processes (anharmonic: Umklapp, 4-phonon) exchange energy with
        the phonon bath, so they rethermalize the mode; elastic processes
        (point-defect/alloy mass-disorder, Rayleigh ~ omega^4) only redirect the
        phonon and conserve its frequency and branch. By default all scattering
        is treated as inelastic; materials with a known elastic component
        override this and derive phonon_relaxation_time from the sum of rates.
        """
        return 1 / self.phonon_relaxation_time(omega), 0.0

    @abstractmethod
    def assign_heat_capacity(self):
        """Calculate heat capacity [J/kg/K] in 3 - 300K range using the polynomial fits"""
        pass

    def group_velocity(self, branch_number, f):
        """
        Group velocity dw/dk [m/s] at ordinary frequency f [Hz] on the given branch,
        as a finite difference between the two tabulated points closest to f.
        Nearest-point search rather than bisection because some branches are not
        monotonic in frequency (e.g. the TA branches of SiC and Graphite).
        """
        f_branch = self.dispersion[:, branch_number + 1]
        diffs = np.abs(f_branch - f)
        nearest = diffs.argmin()
        n = len(f_branch)
        # Interval between the nearest point and whichever of its two neighbors
        # is itself closer to f, so the difference is taken on the side f
        # actually sits on. Using only "nearest - 1" here (as an earlier version
        # did) makes two different phonons whose frequencies straddle the same
        # nearest grid point collapse onto the same interval, giving them the
        # exact same (wrong) velocity instead of two distinct, locally correct
        # ones; comparing both neighbors and picking the closer one avoids that.
        if nearest == 0:
            point_num = 0
        elif nearest == n - 1:
            point_num = nearest - 1
        elif diffs[nearest - 1] <= diffs[nearest + 1]:
            point_num = nearest - 1
        else:
            point_num = nearest
        d_omega = 2 * pi * abs(f_branch[point_num + 1] - f_branch[point_num])
        d_k = abs(self.dispersion[point_num + 1, 0] - self.dispersion[point_num, 0])
        return d_omega / d_k

    def assign_phonon_sampling_tables(self):
        """
        Build the tables for random sampling of a phonon branch and frequency,
        used both for particle emission and for rethermalization at internal
        scattering (Peraud & Hadjiconstantinou, PRB 84, 205331 (2011)).

        The dispersion k-grid is treated as bins; each bin is weighted by the number
        of phonon modes in it (density of states k^2 dk) times the rate at which the
        given process injects energy into those modes. Particles are equal-energy
        bundles representing the deviation from equilibrium, so the per-mode weight
        is built on the mode heat capacity C(w) = dn/dT * hw (not the energy
        spectrum hw*n, which would overweight low frequencies):
        - emission (weight D(w)*C(w)*v(w)): a hot wall emits into mode w at a rate
          proportional to how fast that mode carries energy away, hence the extra
          group-velocity factor (phonon analog of blackbody effusion);
        - scattering (weight D(w)*C(w)/tau_inelastic(w)): only the inelastic
          (anharmonic) channel exchanges energy with the bath, so detailed balance
          is set by its rate alone — modes are absorbed into the bath at rate
          1/tau_inelastic and must be re-emitted at the same rate. Elastic
          (impurity) events conserve the mode and drop out of the balance; a
          phonon redrawn from this table therefore waits, on average, exactly
          tau_inelastic at the drawn frequency before the next redraw (elastic
          events in between don't redraw), keeping the C(w)-shaped population.
        Note that the tables use the actual dispersion, and NOT the Debye
        approximation: the Debye density of states strongly underweights the slow
        zone-edge phonons and overestimates the thermal conductivity.

        Naming convention: attributes ending in _table hold physical quantities of
        each dispersion bin; attributes ending in _probabilities hold cumulative
        probabilities (probability of landing at or below the entry, hence each
        array ends at exactly 1.0), used for inverse sampling: a single uniform
        random number is mapped to a bin via binary search (np.searchsorted).

        Assigned attributes (each a list of 3 per-branch arrays, except the
        branch probabilities, which are single arrays of 3 cumulative values):
        - frequencies_table: bin-midpoint frequency [Hz] of each dispersion bin;
        - group_velocity_table: group velocity dw/dk [m/s] of each bin;
        - emission_frequency_probabilities, scattering_frequency_probabilities:
          cumulative probabilities of the frequency bins within each branch;
        - emission_branch_probabilities, scattering_branch_probabilities:
          cumulative probabilities of the branches (LA, TA1, TA2), for picking
          a branch before picking a frequency.
        """
        self.frequencies_table = []
        self.group_velocity_table = []
        self.emission_frequency_probabilities = []
        self.scattering_frequency_probabilities = []
        emission_branch_weights = []
        scattering_branch_weights = []
        for branch in range(3):
            # Wavevector and frequency grids for this branch
            k_vec = self.dispersion[:, 0]
            f_branch = self.dispersion[:, branch + 1]
            # Bin-midpoint wavevector and bin width in k
            k_mid = (k_vec[1:] + k_vec[:-1]) / 2
            d_k = np.diff(k_vec)
            # Bin-midpoint frequency
            freqs = (f_branch[1:] + f_branch[:-1]) / 2
            # Bin-averaged group velocity dw/dk
            group_velocity = 2 * pi * np.abs(np.diff(f_branch)) / d_k
            # Drop the zero-frequency point at k=0
            valid = freqs > 0
            omegas = 2 * pi * freqs[valid]
            # hw/kT
            x = hbar * omegas / (k_B * self.temp)
            # Mode heat capacity C(w)
            heat_capacity = k_B * x**2 * np.exp(x) / np.expm1(x)**2
            inelastic_rates = np.array([self.phonon_scattering_rates(omega)[0] for omega in omegas])
            # Density of states weight k^2 dk
            dos = k_mid[valid]**2 * d_k[valid]

            self.frequencies_table.append(freqs[valid])
            self.group_velocity_table.append(group_velocity[valid])

            # Cumulative bin weights; the last element is the total weight of
            # the branch, used for the branch probabilities below:
            emission_cumulative = np.cumsum(dos * heat_capacity * group_velocity[valid])
            scattering_cumulative = np.cumsum(dos * heat_capacity * inelastic_rates)
            self.emission_frequency_probabilities.append(emission_cumulative / emission_cumulative[-1])
            self.scattering_frequency_probabilities.append(scattering_cumulative / scattering_cumulative[-1])
            emission_branch_weights.append(emission_cumulative[-1])
            scattering_branch_weights.append(scattering_cumulative[-1])

        self.emission_branch_probabilities = np.cumsum(emission_branch_weights) / np.sum(emission_branch_weights)
        self.scattering_branch_probabilities = np.cumsum(scattering_branch_weights) / np.sum(scattering_branch_weights)

    def assign_dispersion_heat_capacity(self):
        """
        Volumetric heat capacity [J/K/m^3] summed over only the branches present in self.dispersion
        mode C(w) = k*x^2*exp(x)/expm1(x)^2 weighted by the k^2 dk density of states.
        Unlike the experimental fit in assign_heat_capacity, it excludes any physics not in
        the tabulated dispersion (e.g. optical branches), which makes it self-consistent
        with the dispersion-based sampling and the RTA integral.
        """
        k_vec = self.dispersion[:, 0]
        k_mid = (k_vec[1:] + k_vec[:-1]) / 2
        d_k = np.diff(k_vec)
        total = 0.0
        for branch in range(1, self.dispersion.shape[1]):
            freqs = (self.dispersion[1:, branch] + self.dispersion[:-1, branch]) / 2
            valid = freqs > 0
            omegas = 2 * pi * freqs[valid]
            x = hbar * omegas / (k_B * self.temp)
            mode_heat_capacity = k_B * x**2 * np.exp(x) / np.expm1(x)**2
            dos = k_mid[valid]**2 * d_k[valid] / (2 * pi**2)
            total += np.sum(dos * mode_heat_capacity)
        self.dispersion_heat_capacity = total


class Si(Material):
    """
    Physical properties of silicon.
    Dispersion - Ref. Hopkins et al., APL 95, 161902 (2009)
    Relaxation time - impurity and Umklapp coefficients (A, B) fit to bulk
      single-crystal kappa(T), Glassbrenner & Slack, Phys. Rev. 134, A1058
      (1964); see Data/Fitting_BulkSi/ for the fit script and data.
    Heat capacity - Desai P.D. Journal of Physical and Chemical Reference Data 15, 67 (1986)
    Effective mass - H.D. Barber, Effective mass and intrinsic concentration in silicon, Solid-State Electronics, Volume 10, Issue 11 (1967)
    """

    def __init__(self, temp, num_points=1000, fermi_level=None):
        self.name = "Si"
        self.default_speed = 6000   # [m/s] where ??
        self.density = 2330         # [kg/m^3]
        self.vg = 6000              # vitesse de groupe moyenne approx 24/06
        self.temp = temp
        self.assign_electrical_properties(fermi_level)
        self.assign_phonon_dispersion(num_points)
        self.assign_heat_capacity()
        self.assign_dispersion_heat_capacity()
        self.assign_phonon_sampling_tables()

    def assign_phonon_dispersion(self, num_points):
        """Assign phonon dispersion"""

        coefficients_LA = [-9.70e-19, -2.405e-8, 1369.42, 0]
        coefficients_TA = [7.967e-29, 5.674e-19, -7.711e-8, 1081.74, 0]

        self.dispersion = np.zeros((num_points, 4))
        self.dispersion[:, 0] = np.linspace(0, 12e9, num_points)  # Wavevectors
        self.dispersion[:, 1] = np.abs(np.polyval(coefficients_LA, self.dispersion[:, 0]))  # LA branch
        self.dispersion[:, 2] = np.abs(np.polyval(coefficients_TA, self.dispersion[:, 0]))  # TA branch
        self.dispersion[:, 3] = self.dispersion[:, 2]

    def phonon_scattering_rates(self, omega):
        """Umklapp (inelastic) and impurity (elastic) scattering rates [1/s]"""
        deb_temp = 152.0  # Debye temperature normal constant calculated for Si
        rate_impurity = 1.8516e-45 * (omega ** 4)
        rate_umklapp = 1.1026e-19 * (omega ** 2) * self.temp * np.exp(-deb_temp / self.temp)
        return rate_umklapp, rate_impurity

    def phonon_relaxation_time(self, omega):
        """Calculate relaxation time at a given frequency and temperature"""
        inelastic_rate, elastic_rate = self.phonon_scattering_rates(omega)
        return 1 / (inelastic_rate + elastic_rate)

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

    def assign_electrical_properties(self, fermi_level):
        """Assign differents electrical properties to the material."""
        self.effective_electron_dos_mass = 1.18 * electron_mass # [kg] at 300K for pure Si, supposed constant for all temperatures (~1-5% error)
        self.effective_electron_susceptibility_mass = 0.54 * electron_mass
        self.effective_hole_dos_mass = 0.81 * electron_mass
        self.effective_electron_mass = 0.26 * electron_mass
        self.effective_hole_mass = 0.23 * electron_mass # light hole

        if fermi_level:
            self.fermi_level = fermi_level
        else:
            self.fermi_level = -0.037 * electron_volt # [J]


class Vacuum:
    def __init__(self, temp=300):
        self.name = "Vacuum"
        self.density = 0.0
        self.heat_capacity = 0.0
        self.temp = temp
        self.dispersion_table = None

    def get_group_velocity(self, omega, branch_number):
        return 0.0

    def get_dispersion(self, branch_number):
        return None

class SiC(Material):
    """
    Physical properties of silicon carbide
    Dispersion - PRB 50 17054 (1994)
    Relaxation time - Joshi et al, JAP 88, 265 (2000)
    Heat capacity - Collins et al. Journal of Applied Physics 68, 6510 (1990)
    """

    # --- Impurity scattering via Matthiessen's rule (disabled; re-enable when ready) ---
    # 4H-SiC electrical properties:
    # Electron effective masses - Ioffe database / Persson & Lindefelt, PRB 54, 10257 (1996)
    # CT parameters - Roschke & Schwierz, IEEE Trans. Electron Devices 48, 1442 (2001)
    # phonon_limited_electron_mfp = 7e-9  # [m] estimate, needs experimental calibration
    # _ct_mu_max = 947e-4    # [m²/V·s]
    # _ct_mu_min = 40e-4     # [m²/V·s]
    # _ct_n_ref  = 2.0e23    # [m⁻³]  (= 2.0e17 cm⁻³)
    # _ct_alpha  = 0.61
    # ---------------------------------------------------------------------------------

    def __init__(self, temp, num_points=1000, fermi_level=None):
        self.name = "SiC"
        self.density = 3215         # [kg/m^3]
        self.default_speed = 6500   # [m/s] Need to change this probably...
        self.temp = temp
        self.assign_phonon_dispersion(num_points)
        self.assign_heat_capacity()
        self.assign_dispersion_heat_capacity()
        self.assign_phonon_sampling_tables()

    # def assign_electrical_properties(self, fermi_level):
    #     """Assign electrical properties for 4H-SiC"""
    #     self.effective_electron_dos_mass = 1.0 * electron_mass   # [kg] 4H-SiC, Ioffe
    #     self.effective_electron_mass = 0.36 * electron_mass       # [kg] conductivity mass
    #     self.effective_hole_dos_mass = 0.6 * electron_mass
    #     self.effective_hole_mass = 0.6 * electron_mass
    #     self.fermi_level = fermi_level  # None = no CT correction applied

    # def carrier_density(self):
    #     """Electron carrier density from Fermi level using Boltzmann approximation [m⁻³]"""
    #     nc = 2 * (2 * np.pi * self.effective_electron_dos_mass * k_B * self.temp / h_planck**2) ** 1.5
    #     return nc * np.exp(self.fermi_level / (k_B * self.temp))

    # def effective_electron_mfp(self):
    #     """Effective electron MFP combining phonon and impurity scattering via Matthiessen's rule"""
    #     if self.fermi_level is None:
    #         return self.phonon_limited_electron_mfp
    #     n = self.carrier_density()
    #     mu = self._ct_mu_min + (self._ct_mu_max - self._ct_mu_min) / (1 + (n / self._ct_n_ref) ** self._ct_alpha)
    #     return self.phonon_limited_electron_mfp * mu / self._ct_mu_max

    def assign_phonon_dispersion(self, num_points):
        """Assign phonon dispersion"""

        coefficients_LA = [-3.48834e-18, 1.7604452e-08, 1737.36296, 0]
        coefficients_TA = [-2.21696e-19, -3.43668e-08, 1077.98941, 0]

        self.dispersion = np.zeros((num_points, 4))
        self.dispersion[:, 0] = np.linspace(0, 14414281503, num_points)  # Wavevectors
        self.dispersion[:, 1] = np.abs(np.polyval(coefficients_LA, self.dispersion[:, 0]))  # LA branch
        self.dispersion[:, 2] = np.abs(np.polyval(coefficients_TA, self.dispersion[:, 0]))  # TA branch
        self.dispersion[:, 3] = self.dispersion[:, 2]

    def phonon_scattering_rates(self, omega):
        """Umklapp + 4-phonon (inelastic) and impurity (elastic) scattering rates [1/s]"""
        deb_temp = 1200
        rate_impurity = 8.46e-45 * (omega ** 4)
        rate_umklapp = 6.16e-20 * (omega ** 2) * self.temp * np.exp(-deb_temp / self.temp)
        rate_4p = 6.9e-23 * (self.temp ** 2) * (omega ** 2)
        return rate_umklapp + rate_4p, rate_impurity

    def phonon_relaxation_time(self, omega):
        """Calculate relaxation time at a given frequency and temperature including 4 phonon scattering"""
        inelastic_rate, elastic_rate = self.phonon_scattering_rates(omega)
        return 1 / (inelastic_rate + elastic_rate)


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


class Graphite(Material):
    """
    Physical properties of graphite.
    Dispersion - Carbon 91 266-274 (2015)
    Relaxation time - Ref. PRB 87, 115421 (2013)
    Heat capacity - Isaacs, L.L.; Wang, W.Y., Therm. Conduct. 17th, 55-61 (1981)
    """
    dispersion_branch_names = ['LA', 'TA', 'ZA']

    def __init__(self, temp, num_points=1000):
        self.name = "Graphite"
        self.density = 2230            # [kg/m^3]
        self.default_speed = 12900     # [m/s]
        self.temp = temp
        self.assign_phonon_dispersion(num_points)
        self.assign_heat_capacity()
        self.assign_dispersion_heat_capacity()
        self.assign_phonon_sampling_tables()

    def assign_phonon_dispersion(self, num_points):
        """Assign phonon dispersion"""

        coefficients_LA = [-1.24989e-18, -4.11304e-08, 3640.918, 0]
        coefficients_TA = [-1.52298e-18, -4.72535e-08, 2304.367, 0]
        coefficients_ZA = [-3.50255e-18, 1.114371e-07, 106.8250, 0]

        self.dispersion = np.zeros((num_points, 4))
        self.dispersion[:, 0] = np.linspace(0, 14500000000, num_points)  # Wavevectors
        self.dispersion[:, 1] = np.abs(np.polyval(coefficients_LA, self.dispersion[:, 0]))  # LA branch
        self.dispersion[:, 2] = np.abs(np.polyval(coefficients_TA, self.dispersion[:, 0]))  # TA branch
        self.dispersion[:, 3] = np.abs(np.polyval(coefficients_ZA, self.dispersion[:, 0]))  # ZA branch

    def phonon_relaxation_time(self, omega):
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

class SiGe(Material):
    """
    Physical properties of Si0.8Ge0.2 alloy.
    Dispersion - Adapted from Muta et al, J. Alloys Compd. 392, 306 (2005) and Li, J. Phys. Chem. Ref. Data 9, 561 (1980)
    Relaxation time - Umklapp coefficient carried over from the pre-refit Si values
      (Maire et al., Sci. Rep. 7, 41794 (2017) era), impurity/alloy coefficient an
      ad-hoc adaptation of the same; NOT refit against bulk SiGe kappa data, unlike
      Si (see Data/Fitting_BulkSi/). Treat absolute SiGe kappa values with caution.
    Heat capacity - Fit based on Desai P.D., J. Phys. Chem. Ref. Data 13, 1069 (1984) and Wunderlich, Thermophysical Properties of Materials, Springer (2005)
    Effective masses - linearly interpolated between Si and Ge at x=0.2 (Schaffler, 2001)
    """

    def __init__(self, temp, num_points=1000):
        self.name = "SiGe"
        self.default_speed = 3700   # [m/s] – average LA/TA
        self.density = 3008         # [kg/m^3] Si1-xGex: (2.329+3.493x-0.499x**2) g/cm^3, x=0.2, Schaffler (2001)
        self.temp = temp
        self.vg = 3700              # average group velocity approximation 24/06
        self.effective_electron_dos_mass = 1.06 * electron_mass  # [kg] Si0.8Ge0.2: linear interp. Si(1.18) and Ge(0.56)
        self.effective_electron_mass = 0.23 * electron_mass      # [kg] conductivity mass, linear interp. Si(0.26) and Ge(0.12)
        self.effective_hole_dos_mass = 0.71 * electron_mass      # [kg] linear interp. Si(0.81) and Ge(0.29)
        self.effective_hole_mass = 0.19 * electron_mass          # [kg] light hole, linear interp. Si(0.23) and Ge(0.044)
        self.assign_phonon_dispersion(num_points)
        self.assign_heat_capacity()
        self.assign_dispersion_heat_capacity()
        self.assign_phonon_sampling_tables()

    def assign_phonon_dispersion(self, num_points):
        """Assign phonon dispersion"""

        # Coefficients for approximation f(k) from Ge data – need to be changed
        coefficients_LA = [-2.0e-19, -1.0e-8, 1245.0, 0]
        coefficients_TA = [5.0e-29, 4.0e-19, -6.0e-8, 950.0, 0]

        self.dispersion = np.zeros((num_points, 4))
        self.dispersion[:, 0] = np.linspace(0, 12e9, num_points)
        self.dispersion[:, 1] = np.abs(np.polyval(coefficients_LA, self.dispersion[:, 0]))  # LA
        self.dispersion[:, 2] = np.abs(np.polyval(coefficients_TA, self.dispersion[:, 0]))  # TA
        self.dispersion[:, 3] = self.dispersion[:, 2]  # TA2 = TA1 (approx)

    def phonon_scattering_rates(self, omega):
        """Umklapp (inelastic) and impurity/alloy mass-disorder (elastic) scattering rates [1/s].
        The elastic channel dominates in the alloy, so the split matters most here."""
        # Debye Température approx from Ge : ~230 K
        deb_temp = 586.8 #(640 - 266x) K	300 K	Schaffler F. et al.(2001) 4/07
        rate_impurity = 3.5e-45 * (omega ** 4)  # includes alloy mass-disorder; adapted from pre-refit Si values, not fit to SiGe data
        rate_umklapp = 1.1e-19 * (omega ** 2) * self.temp * np.exp(-deb_temp / self.temp)
        return rate_umklapp, rate_impurity

    def phonon_relaxation_time(self, omega):
        """Relaxation time model (impurities + Umklapp scattering)"""
        inelastic_rate, elastic_rate = self.phonon_scattering_rates(omega)
        return 1 / (inelastic_rate + elastic_rate)

    def assign_heat_capacity(self):
        """Empirical polynomial fits for Cp vs T"""

        # Fit from experimental data of Desai + Wunderlich
        if self.temp < 20:
            coeffs = np.array([0.00052, -0.0025, 0.0078])
        elif 20 <= self.temp <= 50:
            coeffs = np.array([-0.001, 0.15, -4.1, 36])
        else:
            coeffs = np.array([-3.2e-6, -4.9e-3, 4.5, -145])

        self.heat_capacity = np.polyval(coeffs, self.temp)


class Diamond(Material):
    """
    Physical properties of diamond
    Dispersion - Ref. PRB 58 12899 (1998)
    """

    def __init__(self, temp, num_points=1000):
        self.name = "Diamond"
        self.density = 3500         # [kg/m^3]
        self.default_speed = 20000  # [m/s]
        self.temp = temp
        self.assign_phonon_dispersion(num_points)

    def assign_phonon_dispersion(self, num_points):
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

    def phonon_relaxation_time(self, omega):
        pass


class AlN(Material):
    """
    Physical properties of AlN.
    Dispersion - Yanagisawa et al, Surface and Interface Analysis 37 133-136 (2005)
    """

    def __init__(self, temp, num_points=1000):
        self.name = "AlN"
        self.density = 3255           # [kg/m^3]
        self.default_speed = 6200     # [m/s]
        self.temp = temp              # [K]
        self.dispersion = np.zeros((num_points, 4))
        self.assign_phonon_dispersion(num_points)

    def assign_phonon_dispersion(self, num_points):
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

    def phonon_relaxation_time(self, omega):
        pass

def get_media_class(material_name: str) -> Material:
    if material_name == "Si":
        return Si
    elif material_name == "SiC":
        return SiC
    elif material_name == "Graphite":
        return Graphite
    else:
        raise Exception(f"Material {material_name} is not supported")
