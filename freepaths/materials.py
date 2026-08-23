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

    def phonon_normal_rate(self, omega):
        """
        Normal (momentum-conserving) three-phonon scattering rate [1/s] at angular
        frequency omega. N-processes keep all three phonons inside the first Brillouin
        zone and conserve total crystal momentum hbar*q, so they are NON-resistive:
        they redistribute momentum among modes rather than destroying it. They are
        therefore deliberately EXCLUDED from phonon_relaxation_time and the kappa_RTA
        integral. Adding 1/tau_N to the Matthiessen sum and randomizing direction at
        the event would make N a second Umklapp, add fake resistance, and underestimate
        kappa at low T. This rate is consumed only by the hydrodynamic tracing path
        (cf.phonon_hydrodynamic), which at an N-event resamples the outgoing mode from
        the drifting Bose-Einstein distribution (conserving the local crystal momentum)
        instead of rethermalizing isotropically. Default: 0 (no N model, pure RTA).
        """
        return 0.0

    def group_velocity(self, branch_number, f):
        """
        Group velocity dw/dk [m/s] at ordinary frequency f [Hz] on the given branch,
        as a finite difference between the two tabulated points closest to f.
        Nearest-point search rather than bisection because some branches are not monotonic in frequency
        """
        f_branch = self.dispersion[:, branch_number + 1]
        diffs = np.abs(f_branch - f)
        nearest = diffs.argmin()
        n = len(f_branch)
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

    def wavevector(self, branch_number, f):
        """Wavevector magnitude k [1/m] at ordinary frequency f [Hz] on the given branch"""
        f_branch = self.dispersion[:, branch_number + 1]
        nearest = np.abs(f_branch - f).argmin()
        return self.dispersion[nearest, 0]

    def mean_inverse_phase_velocity_sq(self):
        """
        Heat-capacity-weighted average of (k/omega)^2 = 1/v_phase^2 over all tabulated
        dispersion modes [s^2/m^2], = sum(DOS*C*(k/omega)^2) / sum(DOS*C). This is the
        material constant that converts the recorded crystal-momentum density into a
        drift velocity in the linearized (small-drift) displaced-Bose-Einstein picture:
        the crystal momentum density is P = chi*u with chi = (T/3) sum(DOS*C*(k/omega)^2),
        and the deviational energy density is e = C_v*T with C_v = sum(DOS*C)
        (dispersion_heat_capacity), so
            u = 3 * (P / e) / <(k/omega)^2>_C ,
        with T, volume and sample count cancelling. Dominated by the slow, large-k/omega
        modes (the quadratic ZA branch), consistent with the flexural branch carrying the
        phonon hydrodynamics. Cached after the first call.
        """
        if getattr(self, "_mean_inv_vp2", None) is not None:
            return self._mean_inv_vp2
        k_vec = self.dispersion[:, 0]
        k_mid = (k_vec[1:] + k_vec[:-1]) / 2
        d_k = np.diff(k_vec)
        numerator = 0.0
        denominator = 0.0
        for branch in range(1, self.dispersion.shape[1]):
            freqs = (self.dispersion[1:, branch] + self.dispersion[:-1, branch]) / 2
            valid = freqs > 0
            omegas = 2 * pi * freqs[valid]
            x = hbar * omegas / (k_B * self.temp)
            mode_heat_capacity = k_B * x**2 * np.exp(x) / np.expm1(x)**2
            dos = k_mid[valid]**2 * d_k[valid]
            inv_vp2 = (k_mid[valid] / omegas)**2
            numerator += np.sum(dos * mode_heat_capacity * inv_vp2)
            denominator += np.sum(dos * mode_heat_capacity)
        self._mean_inv_vp2 = numerator / denominator
        return self._mean_inv_vp2

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
        Unlike an experimental heat-capacity fit, it excludes any physics not in the tabulated
        dispersion (e.g. optical branches), which makes it self-consistent with the
        dispersion-based sampling and the RTA integral. This is the heat capacity used for the
        temperature-profile conversion (maps.py).
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

    def ballistic_conductance(self):
        """
        Landauer ballistic thermal conductance per unit cross-sectional area [W/m^2/K]:
            G = (1/4) sum_branches int (d^3k / (2 pi)^3) v(w) C(w),
        where C(w) = hbar*w * df_eq/dT is the mode heat capacity and the 1/4 is the
        isotropic forward-flux projection <v_x * Theta(v_x)> = v/4.
        Huang et al., Nat. Commun. 14, 2044 (2023), Eq. 1
        """
        k_vec = self.dispersion[:, 0]
        k_mid = (k_vec[1:] + k_vec[:-1]) / 2
        d_k = np.diff(k_vec)
        total = 0.0
        for branch in range(1, self.dispersion.shape[1]):
            freqs = (self.dispersion[1:, branch] + self.dispersion[:-1, branch]) / 2
            group_velocity = 2 * pi * np.abs(np.diff(self.dispersion[:, branch])) / d_k
            valid = freqs > 0
            omegas = 2 * pi * freqs[valid]
            x = hbar * omegas / (k_B * self.temp)
            mode_heat_capacity = k_B * x**2 * np.exp(x) / np.expm1(x)**2
            dos = k_mid[valid]**2 * d_k[valid] / (2 * pi**2)
            total += np.sum(dos * mode_heat_capacity * group_velocity[valid]) / 4.0
        return total


class Si(Material):
    """
    Physical properties of silicon.
    Dispersion - Ref. Hopkins et al., APL 95, 161902 (2009)
    Impurity and Umklapp coefficients (A, B) fit to bulk kappa(T), Glassbrenner & Slack, Phys. Rev. 134, A1058 (1964)
    Heat capacity - Desai P.D. Journal of Physical and Chemical Reference Data 15, 67 (1986)
    Effective mass - H.D. Barber, Effective mass and intrinsic concentration in silicon, Solid-State Electronics, Volume 10, Issue 11 (1967)
    """

    def __init__(self, temp, num_points=1000, fermi_level=None):
        self.name = "Si"
        self.density = 2330         # [kg/m^3]
        self.vg = 6000              # vitesse de groupe moyenne approx 24/06
        self.temp = temp
        self.assign_electrical_properties(fermi_level)
        self.assign_phonon_dispersion(num_points)
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

    def assign_electrical_properties(self, fermi_level):
        """Assign differents electrical properties to the material."""
        self.effective_electron_dos_mass = 1.18 * electron_mass # [kg] at 300K for pure Si
        self.effective_electron_susceptibility_mass = 0.54 * electron_mass
        self.effective_hole_dos_mass = 0.81 * electron_mass
        self.effective_electron_mass = 0.26 * electron_mass
        self.effective_hole_mass = 0.23 * electron_mass # light hole
        self.dielectric_constant = 11.7  # static relative permittivity of Si

        if fermi_level:
            self.fermi_level = fermi_level
        else:
            self.fermi_level = -0.037 * electron_volt # [J]


class Vacuum:
    def __init__(self, temp=300):
        self.name = "Vacuum"
        self.density = 0.0
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
        self.temp = temp
        self.assign_phonon_dispersion(num_points)
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



class Graphite(Material):
    """
    Physical properties of graphite.
    Dispersion - Carbon 91 266-274 (2015)
    Relaxation time - Ref. PRB 87, 115421 (2013)
    Heat capacity - Isaacs, L.L.; Wang, W.Y., Therm. Conduct. 17th, 55-61 (1981)
    """
    dispersion_branch_names = ['LA', 'TA', 'ZA']

    def __init__(self, temp, num_points=1000, isotope_c13_concentration=0.0):
        self.name = "Graphite"
        self.density = 2230            # [kg/m^3]
        self.temp = temp
        self.assign_phonon_dispersion(num_points)
        # Isotope table before the sampling tables, which query phonon_scattering_rates:
        self.assign_isotope_scattering(isotope_c13_concentration)
        self.assign_dispersion_heat_capacity()
        self.assign_phonon_sampling_tables()

    def assign_phonon_dispersion(self, num_points):
        """Assign phonon dispersion"""

        coefficients_LA = [-1.24989e-18, -4.11304e-08, 3640.918, 0]
        coefficients_TA = [-1.52298e-18, -4.72535e-08, 2304.367, 0]

        self.dispersion = np.zeros((num_points, 4))
        self.dispersion[:, 0] = np.linspace(0, 14500000000, num_points)  # Wavevectors
        self.dispersion[:, 1] = np.abs(np.polyval(coefficients_LA, self.dispersion[:, 0]))  # LA branch
        self.dispersion[:, 2] = np.abs(np.polyval(coefficients_TA, self.dispersion[:, 0]))  # TA branch

        # ZA (flexural) branch: quadratic near Gamma, omega = b_ZA * k^2 (Nihira & Iwata
        # semicontinuum model via Alofi & Srivastava, PRB 87, 115421 (2013), Eqs. 5/13;
        # bending parameter b = 3.13e-3 cm^2/s = 3.13e-7 m^2/s). 
        b_ZA = 3.13e-7  # [m^2/s]
        self.dispersion[:, 3] = b_ZA * self.dispersion[:, 0] ** 2 / (2 * pi)  # ZA branch [Hz]

    # Three-phonon Umklapp rate, Klemens/Slack high-T form
    #   1/tau_U = B_U * omega^2 * T * exp(-deb_temp/(alpha*T))
    _B_N = 2.12e-25      # Normal-process prefactor (momentum-conserving) [Alofi, unused; see phonon_normal_rate]
    _B_U = 4.003e-21     # Umklapp prefactor (resistive), Slack omega^2*T*exp form, fit to bulk kappa via Callaway
    _deb_temp = 100.0    # effective Umklapp activation temperature theta* [K] (with alpha=1); NOT the thermodynamic Debye temp
    _alpha = 1.0         # (theta* = deb_temp/alpha)

    # Normal-process prefactor [Hz/K^3], calibrated so the high-frequency plateau reaches the
    # first-principles graphite N-rate ~1e10 Hz at 100 K (Guo et al., PRB 104, 075450 (2021),
    # Fig. 1a). 
    _C_N = 1.0e4
    # Normal-process frequency shape g(w) = (w/wc)^p / (1 + (w/wc)^p): a saturating rise that
    # suppresses the low-frequency N and preserves the high-f plateau. Fit to Guo Fig. 1a (-m)
    # Normal rate at 100 K (which falls to ~1e9 at low f, 10x below the plateau); it also brings
    # l_N(60 K) to ~4-11 um (Huang SI Fig. 6, ~2-5 um) instead of the old flat-rate ~1 um.
    _N_shape_wc = 2 * pi * 2.93e12   # crossover angular frequency [rad/s]
    _N_shape_p = 1.38                # crossover sharpness

    # 13-C isotope masses [amu] for the mass-variance parameter (12-C defines the amu):
    _M_C12 = 12.0
    _M_C13 = 13.003355

    def assign_isotope_scattering(self, concentration):
        """
        Precompute the Tamura point-defect (isotope mass-disorder) scattering-rate table
        1/tau_iso(w) = (pi/6) g^2 w^2 D(w) from the tabulated dispersion, where
        g^2 = sum_i c_i (1 - m_i/<m>)^2 is the mass variance (13-C at fractional abundance
        `concentration` in carbon; natural = 0.0107 -> g^2 ~ 7.4e-5) and D(w) is the total
        phonon DOS per unit cell, built by histogramming the k^2 dk mode weight over the
        three branches and normalized so each branch integrates to one mode/cell (the
        Tamura final-state sum rule). That normalization makes the rate independent of the
        fitted BZ k-extent and, in the Debye limit D -> 3 V0 w^2 / (2 pi^2 v^3), recovers
        the Klemens form 1/tau = V0 g^2 w^4 / (4 pi v^3). 
        Tamura, PRB 27, 858 (1983); Klemens, Proc. Phys. Soc. A 68, 1113 (1955).
        """
        c13 = concentration
        c12 = 1.0 - c13
        m_avg = c12 * self._M_C12 + c13 * self._M_C13
        self._isotope_g2 = c12 * (1 - self._M_C12 / m_avg) ** 2 + c13 * (1 - self._M_C13 / m_avg) ** 2
        if self._isotope_g2 == 0.0:
            self._isotope_f_grid = None
            self._isotope_rate_grid = None
            return
        # Total DOS per unit cell, summed over branches, on a uniform frequency histogram.
        # Each k-bin holds 3 k^2 dk / k_max^3 modes/cell (so each branch integrates to 1);
        # dividing the binned mode count by the bin's angular-frequency width gives D(w) [s].
        k_vec = self.dispersion[:, 0]
        k_max = k_vec[-1]
        k_mid = (k_vec[1:] + k_vec[:-1]) / 2
        d_k = np.diff(k_vec)
        mode_weight = 3.0 * k_mid ** 2 * d_k / k_max ** 3
        f_max = self.dispersion[:, 1:4].max()
        n_bins = 400
        f_edges = np.linspace(0, f_max, n_bins + 1)
        modes_per_cell = np.zeros(n_bins)
        for branch in range(1, 4):
            f_mid = (self.dispersion[1:, branch] + self.dispersion[:-1, branch]) / 2
            modes_per_cell += np.histogram(f_mid, bins=f_edges, weights=mode_weight)[0]
        d_omega = 2 * pi * (f_edges[1] - f_edges[0])
        dos_per_cell = modes_per_cell / d_omega                       # D(w) [s], total integral = 3
        omega_grid = 2 * pi * (f_edges[1:] + f_edges[:-1]) / 2
        self._isotope_f_grid = omega_grid / (2 * pi)
        self._isotope_rate_grid = (pi / 6) * self._isotope_g2 * omega_grid ** 2 * dos_per_cell

    def phonon_isotope_rate(self, omega):
        """Elastic isotope (mass-disorder) scattering rate [1/s] at angular frequency omega,
        interpolated from the precomputed Tamura table; 0 for an isotopically pure crystal."""
        if self._isotope_rate_grid is None:
            return 0.0
        return float(np.interp(omega / (2 * pi), self._isotope_f_grid, self._isotope_rate_grid))

    def phonon_scattering_rates(self, omega):
        """Return (inelastic, elastic) scattering rates [1/s] """
        rate_umklapp = self._B_U * (omega ** 2) * self.temp * np.exp(-self._deb_temp / (self._alpha * self.temp))
        return rate_umklapp, self.phonon_isotope_rate(omega)

    def phonon_relaxation_time(self, omega):
        """Calculate relaxation time at a given frequency and temperature"""
        inelastic_rate, elastic_rate = self.phonon_scattering_rates(omega)
        return 1 / (inelastic_rate + elastic_rate)

    def phonon_normal_rate(self, omega):
        """
        Normal (momentum-conserving) three-phonon rate [1/s]: Callaway T^3 scaling with a
        saturating frequency shape, 1/tau_N = C_N * T^3 * (w/wc)^p / (1 + (w/wc)^p).
        The shape -> 1 at high frequency (the ~1e10 Hz plateau at 100 K, Guo et al., PRB 104, 075450
        (2021), Fig. 1a) and suppresses the low-frequency N, where the ab-initio Normal rate
        falls to ~1e9 (10x below the plateau). This brings the Normal mean free path at 60 K to
        ~4-11 um (Huang et al., Nat. Commun. 14, 2044 (2023), SI Fig. 6, ~2-5 um), versus the
        old frequency-INDEPENDENT 1/tau_N = C_N*T^3 which gave l_N ~ 1 um at low f (too
        collisional) and thus an over-driven drift / over-inflated hydrodynamic hump. N is
        momentum-CONSERVING, so this does NOT affect kappa_RTA; it only softens the drift.
        """
        shape = (omega / self._N_shape_wc) ** self._N_shape_p
        return self._C_N * (self.temp ** 3) * shape / (1.0 + shape)



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
        self.density = 3008         # [kg/m^3] Si1-xGex: (2.329+3.493x-0.499x**2) g/cm^3, x=0.2, Schaffler (2001)
        self.temp = temp
        self.vg = 3700              # average group velocity approximation 24/06
        self.effective_electron_dos_mass = 1.06 * electron_mass  # [kg] Si0.8Ge0.2: linear interp. Si(1.18) and Ge(0.56)
        self.effective_electron_mass = 0.23 * electron_mass      # [kg] conductivity mass, linear interp. Si(0.26) and Ge(0.12)
        self.effective_hole_dos_mass = 0.71 * electron_mass      # [kg] linear interp. Si(0.81) and Ge(0.29)
        self.effective_hole_mass = 0.19 * electron_mass          # [kg] light hole, linear interp. Si(0.23) and Ge(0.044)
        self.assign_phonon_dispersion(num_points)
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



class Diamond(Material):
    """
    Physical properties of diamond
    Dispersion - Ref. PRB 58 12899 (1998)
    """

    def __init__(self, temp, num_points=1000):
        self.name = "Diamond"
        self.density = 3500         # [kg/m^3]
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
