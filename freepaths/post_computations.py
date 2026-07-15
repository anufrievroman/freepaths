import numpy as np
from scipy.constants import k, electron_volt, elementary_charge, hbar
from scipy.integrate import simpson

from freepaths.data import GeneralData
from freepaths.materials import get_media_class
from freepaths.config import cf

class ElectronPostComputation:
    """Handle computations for electrons"""

    def __init__(self, general_data: GeneralData):
        """General_data must have collected data from workers before initialization"""
        self.electron_data = np.column_stack((general_data.initial_energies, general_data.travel_times))
        self.energies = self.electron_data[:, 0]
        self.travel_times = self.electron_data[:, 1]
        self.cold_side_mask = None
        self.energies_unique = None
        self.energies_inv_index = None
        self.mean_travel_times = None
        self.mc_tdf = None
        self.mc_conductivity = None
        self.mc_seebeck = None
        self.mc_power_factor = None
        self.mc_thermal_conductivity = None
        self.mc_zt = None
        self.mapping_constant = None

        material_class = get_media_class(cf.media)
        self.material = material_class(cf.temp, fermi_level=cf.media_fermi_level)
        self.effective_dos_mass = None
        self.effective_conduction_mass = None

        self.fermi_levels = None
        self.fermi_dist = None
        self.fermi_prime_dist = None
        self.density_of_states = None
        self.eta = None
        self.bte_tdf = None
        self.bte_conductivity = None
        self.bte_seebeck = None
        self.bte_power_factor = None
        self.bte_thermal_conductivity = None
        self.mean_mapping_constant = None

    def compute_unique_values(self):
        """Compute unique energies values that crossed the cold side"""
        self.cold_side_mask = self.electron_data[:, 1] != 0 # Keep only electrons that crossed the cold side
        energies = self.electron_data[self.cold_side_mask, 0]

        # Get indices for each electron with the same energy
        self.energies_unique, self.energies_inv_index = np.unique(energies, return_inverse=True)

    def compute_physical_functions(self):
        """Compute physical functions used for others calculations"""
        # Resolution of the Fermi level sweep follows the electron energy step:
        n_points = round((cf.fermi_level_upper_bound - cf.fermi_level_lower_bound) / cf.energy_step)
        self.fermi_levels = np.linspace(cf.fermi_level_lower_bound * electron_volt,
                                        cf.fermi_level_upper_bound * electron_volt, n_points)

        # η = (E - Ef) / kT: dimensionless energy measured from the Fermi level, used in transport integrals
        self.eta = np.subtract.outer(self.energies_unique, self.fermi_levels) / (k*cf.temp) # shape = (energies_unique, fermi_levels)

        # Fermi-Dirac distribution: f(E) = 1 / (1 + exp(η))
        self.fermi_dist = 1.0/(1.0 + np.exp(self.eta))

        # Energy derivative of f: ∂f/∂E = -f(1-f) / kT (acts as a window around Ef)
        self.fermi_prime_dist = (- self.fermi_dist * (1-self.fermi_dist) / (k*cf.temp))

        if cf.is_carrier_electron:
            self.effective_dos_mass = self.material.effective_electron_dos_mass
            self.effective_conduction_mass = self.material.effective_electron_mass
        else:
            self.effective_dos_mass = self.material.effective_hole_dos_mass
            self.effective_conduction_mass = self.material.effective_hole_mass

        # 3D parabolic band DoS: g(E) = (1/2π²) × (2m*/ℏ²)^(3/2) × √E
        self.density_of_states = (1/(2*(np.pi**2))) * ((2*self.effective_dos_mass/(hbar**2))**(3/2)) * np.sqrt(self.energies_unique)

    def compute_mean_travel_times(self):
        """Compute mean travel times with respect to initial energy (Priyadarshi et al. 2023, Eq. 1)"""
        travel_times = self.electron_data[self.cold_side_mask, 1]

        # Sum together travel_times for the same energy
        sums = np.bincount(self.energies_inv_index, weights=travel_times)
        # Get the number of electron for a given energy
        counts = np.bincount(self.energies_inv_index)
        means = sums/counts
        self.mean_travel_times = np.column_stack((self.energies_unique, means))

    def compute_bte_tdf(self):
        """Compute analytical BTE TDF for pristine material with constant MFP"""
        mfp = cf.electron_mfp
        # From BTE: Ξ_BTE = (1/3)τv²g; the 1/3 is the 3D angular average ⟨v_x²/v²⟩; with τ = mfp/v: Ξ_BTE = (mfp/3)·v(E)·g(E)
        self.bte_tdf = np.column_stack((self.energies_unique, (mfp/3) * (2*self.energies_unique/(self.effective_conduction_mass))**0.5 * self.density_of_states))

    def compute_mapping_constant(self):
        # C = σ_BTE(Ef) / σ_MC_raw(Ef): calibrates the MC flux to physical units at the material's Fermi level
        mc_raw_tdf = (1.0 / self.mean_travel_times[:,1]) * self.density_of_states

        integrand = self.bte_tdf[:,1][:,None] * (-self.fermi_prime_dist)
        bte_conductivity = elementary_charge**2 * simpson(integrand, x=self.energies_unique, axis=0)

        integrand = mc_raw_tdf[:,None] * (-self.fermi_prime_dist)
        mc_conductivity = elementary_charge**2 * simpson(integrand, x=self.energies_unique, axis=0)

        C = bte_conductivity / mc_conductivity

        self.mapping_constant = np.column_stack((self.fermi_levels, C))
        self.mean_mapping_constant = np.mean(C)

    def compute_mc_tdf(self):
        """Compute the transport distribution function with respect to initial energy (Priyadarshi et al. 2023, Eq. 3)"""
        if cf.mean_mapping_constant is not None:
            area_constant = cf.mean_mapping_constant # in m^2
        else:
            area_constant = self.mean_mapping_constant
        # Ξ(E) = C × F(E) × g(E), where F(E) = 1/⟨ToF(E)⟩ is the electron flux per simulated particle
        self.mc_tdf = np.column_stack((self.energies_unique, area_constant * (1/self.mean_travel_times[:,1]) * self.density_of_states))

        self.mc_tdf_vals = self.mc_tdf[:, 1][:, None] # Shape = (energies_unique, 1)

    def compute_conductivity(self):
        """Compute electron conductivity with respect to fermi-level (Priyadarshi et al. 2023, Eq. 4)"""
        # σ = e² ∫ Ξ(E) × (-∂f/∂E) dE
        integrand = self.bte_tdf[:,1][:,None] * (-self.fermi_prime_dist)
        self.bte_conductivity = elementary_charge**2 * simpson(integrand, x=self.energies_unique, axis=0)
        self.bte_conductivity = np.column_stack((self.fermi_levels, self.bte_conductivity))

        integrand = self.mc_tdf_vals * (-self.fermi_prime_dist)
        self.mc_conductivity = elementary_charge**2 * simpson(integrand, x=self.energies_unique, axis=0)
        self.mc_conductivity = np.column_stack((self.fermi_levels, self.mc_conductivity))

    def compute_seebeck(self):
        """Compute Seebeck coefficient with respect to fermi-level (Priyadarshi et al. 2023, Eq. 5)"""
        # S = ∓(e·kB/σ) ∫ Ξ(E) × (-∂f/∂E) × η dE; the sign follows the carrier charge:
        # negative for electrons (n-type), positive for holes (p-type)
        sign = -1.0 if cf.is_carrier_electron else 1.0
        integrand = self.bte_tdf[:,1][:,None] * (-self.fermi_prime_dist) * self.eta
        self.bte_seebeck = sign * (elementary_charge * k / self.bte_conductivity[:,1]) * simpson(integrand, x=self.energies_unique, axis=0)
        self.bte_seebeck = np.column_stack((self.fermi_levels, self.bte_seebeck))

        integrand = self.mc_tdf_vals * (-self.fermi_prime_dist) * self.eta
        self.mc_seebeck = sign * (elementary_charge * k / self.mc_conductivity[:,1]) * simpson(integrand, x=self.energies_unique, axis=0)
        self.mc_seebeck = np.column_stack((self.fermi_levels, self.mc_seebeck))

    def compute_power_factor(self):
        """Compute the thermoelectric Power Factor (Priyadarshi et al. 2023, Eq. 6): PF = σS²"""
        self.bte_power_factor = np.column_stack((self.fermi_levels, self.bte_conductivity[:,1] * (self.bte_seebeck[:,1]**2)))
        self.mc_power_factor = np.column_stack((self.fermi_levels, self.mc_conductivity[:,1] * (self.mc_seebeck[:,1]**2)))

    def compute_thermal_conductivity(self):
        """Compute the electronic thermal conductivity with respect to fermi-level (Priyadarshi et al. 2023, Eq. 7)"""
        # κe = (1/T) ∫ Ξ(E) × (-∂f/∂E) × (E - Ef)² dE - S²σT; (E - Ef) = η·kT
        integrand = self.bte_tdf[:,1][:,None] * (-self.fermi_prime_dist) * (self.eta * k*cf.temp)**2
        self.bte_thermal_conductivity = (1/cf.temp) * simpson(integrand, x=self.energies_unique, axis=0) - self.bte_power_factor[:,1]*cf.temp
        self.bte_thermal_conductivity = np.column_stack((self.fermi_levels, self.bte_thermal_conductivity))

        integrand = self.mc_tdf_vals * (-self.fermi_prime_dist) * (self.eta * k*cf.temp)**2
        self.mc_thermal_conductivity = (1/cf.temp) * simpson(integrand, x=self.energies_unique, axis=0) - self.mc_power_factor[:,1]*cf.temp
        self.mc_thermal_conductivity = np.column_stack((self.fermi_levels, self.mc_thermal_conductivity))

    def compute_zt(self, kappa_ph):
        """Compute ZT vs Fermi level: ZT = S²σT / (κₑ + κ_ph)"""
        kappa_total = self.mc_thermal_conductivity[:, 1] + kappa_ph
        zt = self.mc_power_factor[:, 1] * cf.temp / kappa_total
        self.mc_zt = np.column_stack((self.fermi_levels, zt))

    def compute_relaxation_time(self):
        # τ(E) = mfp / v(E): relaxation time from constant MFP and parabolic-band group velocity
        self.speeds = ((2 * self.energies_unique) / (self.effective_conduction_mass)) ** 0.5
        self.relaxation_times = cf.electron_mfp / self.speeds
        self.relaxation_times = np.column_stack((self.energies_unique, self.relaxation_times))

    def write_into_file(self):
        np.savetxt("Data/Mean travel time vs energy.csv", self.mean_travel_times, fmt='%2.4e', header="Energy [J], Travel time [s]", encoding='utf-8', delimiter=',')
        tdf_combined = np.column_stack((self.mc_tdf, self.bte_tdf[:, 1]))
        np.savetxt("Data/Transport distribution function.csv", tdf_combined, fmt='%2.4e', header="Energy [J], MC TDF [s^-1 J^-1 m^-1], BTE TDF [s^-1 J^-1 m^-1]", encoding='utf-8', delimiter=',')
        np.savetxt("Data/Electron conductivity.csv", np.column_stack((self.mc_conductivity, self.bte_conductivity[:,1])), fmt='%2.4e', header="Fermi-level [J], MC Conductivity [S/m], BTE Conductivity [S/m]", encoding='utf-8', delimiter=',')
        np.savetxt("Data/Seebeck coefficient.csv", np.column_stack((self.mc_seebeck, self.bte_seebeck[:,1])), fmt='%2.4e', header="Fermi-level [J], MC Seebeck coefficient [V/K], BTE Seebeck coefficient [V/K]", encoding='utf-8', delimiter=',')
        np.savetxt("Data/Power factor.csv", np.column_stack((self.mc_power_factor, self.bte_power_factor[:,1])), fmt='%2.4e', header="Fermi-level [J], MC Power factor [W/m/K^2], BTE Power factor [W/m/K^2]", encoding='utf-8', delimiter=',')
        np.savetxt("Data/Mapping constant.csv", self.mapping_constant, fmt='%2.4e', header="Fermi-level [J], Mapping constant [m^2]", encoding='utf-8', delimiter=',')
        np.savetxt("Data/Thermal conductivity el.csv", np.column_stack((self.mc_thermal_conductivity, self.bte_thermal_conductivity[:,1])), fmt='%2.4e', header="Fermi-level [J], MC Electron thermal conductivity [W/m/K], BTE Electron thermal conductivity [W/m/K]", encoding='utf-8', delimiter=',')
        np.savetxt("Data/Scattering time vs energy.csv", self.relaxation_times, fmt='%2.4e', header="Energy [J], Scattering time [s]", encoding='utf-8', delimiter=',')
        if self.mc_zt is not None:
            np.savetxt("Data/ZT.csv", self.mc_zt, fmt='%2.4e', header="Fermi-level [J], ZT", encoding='utf-8', delimiter=',')

    def compute(self):
        self.compute_unique_values()
        self.compute_physical_functions()
        self.compute_mean_travel_times()
        self.compute_bte_tdf()
        self.compute_mapping_constant()
        self.compute_mc_tdf()
        self.compute_conductivity()
        self.compute_seebeck()
        self.compute_power_factor()
        self.compute_thermal_conductivity()
        self.compute_relaxation_time()
