import numpy as np
from scipy.constants import k, electron_volt, elementary_charge, hbar
from scipy.integrate import simpson

from freepaths.data import GeneralData
from freepaths.materials import get_media_class
from freepaths.config import cf

class ElectronPostComputation:
    """
    Handle computations for electrons
    """
    
    def __init__(self, general_data: GeneralData):
        """
        general_data must have collected data from workers before initialization
        """
        self.electron_data = np.column_stack((general_data.initial_energies, general_data.travel_times))
        self.energies = self.electron_data[:, 0]
        self.travel_times = self.electron_data[:, 1]
        self.cold_side_mask = None
        self.energies_unique = None
        self.energies_inv_index = None
        self.mean_travel_times = None
        self.transport_function = None
        self.electron_conductivity = None
        self.seebeck_coeff = None
        self.power_factor = None
        self.electron_thermal_conductivity = None
        self.figure_of_merit = None
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
        self.true_tdf = None
        self.true_conductivity = None
        self.true_seebeck = None
        self.true_power_factor = None
        self.true_thermal_conductivity = None
        self.true_figure_of_merit = None
        self.mean_mapping_constant = None
        
    def compute_unique_values(self):
        """Compute unique energies values that crossed the cold side"""
        self.cold_side_mask = self.electron_data[:, 1] != 0 # Keep only electrons that crossed the cold side
        energies = self.electron_data[self.cold_side_mask, 0]
        
        # Get indices for each electron with the same energy
        self.energies_unique, self.energies_inv_index = np.unique(energies, return_inverse=True)
    
    def compute_physical_functions(self):
        """Compute physical functions used for others calculations"""
        self.fermi_levels = np.linspace(self.material.fermi_level - 0.1*electron_volt, self.material.fermi_level + 0.1*electron_volt, 200)
        # Suppose energy null at minimum of conduction band
        self.eta = np.subtract.outer(self.energies_unique, self.fermi_levels) / (k*cf.temp) # shape = (energies_unique, fermi_levels)
        self.fermi_dist = 1.0/(1.0 + np.exp(self.eta))
        self.fermi_prime_dist = (- self.fermi_dist * (1-self.fermi_dist) / (k*cf.temp))
        
        if cf.is_carrier_electron:
            self.effective_dos_mass = self.material.effective_electron_dos_mass
            self.effective_conduction_mass = self.material.effective_electron_mass
        else:
            self.effective_dos_mass = self.material.effective_hole_dos_mass
            self.effective_conduction_mass = self.material.effective_hole_mass
        
        self.density_of_states = (1/(2*(np.pi**2))) * ((2*self.effective_dos_mass/(hbar**2))**(3/2)) * np.sqrt(self.energies_unique)
    
    def compute_mean_travel_times(self):
        """Compute mean travel times with respect to initial energy"""
        travel_times = self.electron_data[self.cold_side_mask, 1]
        
        # Sum together travel_times for the same energy
        sums = np.bincount(self.energies_inv_index, weights=travel_times)
        # Get the number of electron for a given energy
        counts = np.bincount(self.energies_inv_index)
        means = sums/counts
        self.mean_travel_times = np.column_stack((self.energies_unique, means))
    
    def compute_transport_function(self):
        """Compute the transport distribution function with respect to initial energy"""
        if cf.mean_mapping_constant is not None:
            area_constant = cf.mean_mapping_constant # in m^2
        else:
            area_constant = self.mean_mapping_constant
        self.transport_function = np.column_stack((self.energies_unique, area_constant * (1/self.mean_travel_times[:,1])**(0.25) * self.density_of_states))
        
        self.Sigma = self.transport_function[:, 1][:, None] # Shape = (energies_unique, 1)
    
    def compute_electron_conductivity(self):
        """Compute electron conductivity with respect to fermi-level"""
        integrand = self.true_tdf[:,1][:,None] * (-self.fermi_prime_dist)
        self.true_conductivity = elementary_charge**2 * simpson(integrand, x=self.energies_unique, axis=0)
        self.true_conductivity = np.column_stack((self.fermi_levels, self.true_conductivity))
            
        integrand = self.Sigma * (-self.fermi_prime_dist)
        self.electron_conductivity = elementary_charge**2 * simpson(integrand, x=self.energies_unique, axis=0)
        self.electron_conductivity = np.column_stack((self.fermi_levels, self.electron_conductivity))
    
    def compute_seebeck(self):
        """Compute Seebeck coefficient with respect to fermi-level"""
        integrand = self.true_tdf[:,1][:,None] * (-self.fermi_prime_dist) * self.eta
        self.true_seebeck = (elementary_charge * k / self.true_conductivity[:,1]) * simpson(integrand, x=self.energies_unique, axis=0)
        self.true_seebeck = np.column_stack((self.fermi_levels, self.true_seebeck))
        
        integrand = self.Sigma * (-self.fermi_prime_dist) * self.eta
        self.seebeck_coeff = (elementary_charge * k / self.electron_conductivity[:,1]) * simpson(integrand, x=self.energies_unique, axis=0)
        self.seebeck_coeff = np.column_stack((self.fermi_levels, self.seebeck_coeff))
    
    def compute_power_factor(self):
        """Compute the thermoelectric Power Factor"""
        self.true_power_factor = np.column_stack((self.fermi_levels, self.true_conductivity[:,1] * (self.true_seebeck[:,1]**2)))
        self.power_factor = np.column_stack((self.fermi_levels, self.electron_conductivity[:,1] * (self.seebeck_coeff[:,1]**2)))
    
    def compute_electronic_thermal_conductivity(self):
        """Compute the electronic thermal conductivity with respect to fermi-level"""
        integrand = self.true_tdf[:,1][:,None] * (-self.fermi_prime_dist) * (self.eta * k*cf.temp)**2
        self.true_thermal_conductivity = (1/cf.temp) * simpson(integrand, x=self.energies_unique, axis=0) - self.true_power_factor[:,1]*cf.temp
        self.true_thermal_conductivity = np.column_stack((self.fermi_levels, self.true_thermal_conductivity))      
        
        
        integrand = self.Sigma * (-self.fermi_prime_dist) * (self.eta * k*cf.temp)**2
        self.electron_thermal_conductivity = (1/cf.temp) * simpson(integrand, x=self.energies_unique, axis=0) - self.power_factor[:,1]*cf.temp
        self.electron_thermal_conductivity = np.column_stack((self.fermi_levels, self.electron_thermal_conductivity))
    
    def compute_figure_of_merit(self):
        """Compute the figure of merit with respect to fermi-level"""
        self.true_figure_of_merit = np.column_stack((self.fermi_levels, self.true_power_factor[:,1] * cf.temp / self.true_thermal_conductivity[:,1]))
        self.figure_of_merit = np.column_stack((self.fermi_levels, self.power_factor[:,1] * cf.temp / self.electron_thermal_conductivity[:,1]))
    
    def compute_true_tdf(self):
        """Compute true tdf for a pristine material, not correct in most other cases"""
        mfp = cf.electron_mfp
        self.true_tdf = np.column_stack((self.energies_unique, mfp * (2*self.energies_unique/(self.effective_conduction_mass))**0.5 * self.density_of_states))
      

    def compute_mapping_constant(self):
        raw_tdf  = (1.0 / self.mean_travel_times[:,1]) ** (0.25) * self.density_of_states

        integrand = self.true_tdf[:,1][:,None] * (-self.fermi_prime_dist)
        true_conductivity = elementary_charge**2 * simpson(integrand, x=self.energies_unique, axis=0)
        
        integrand = raw_tdf[:,None] * (-self.fermi_prime_dist)
        electron_conductivity = elementary_charge**2 * simpson(integrand, x=self.energies_unique, axis=0)
        
        C = true_conductivity / electron_conductivity
        C_fit = np.ones_like(C) * C.mean()
        
        self.mapping_constant = np.column_stack((self.fermi_levels, C, C_fit))
        self.mean_mapping_constant = self.mapping_constant[:,1].mean()
        self.mean_mapping_constant = np.interp(self.material.fermi_level, self.mapping_constant[:,0], self.mapping_constant[:,1])

    def compute_scattering_rate(self):
        self.speeds = ((2 * self.energies_unique) / (self.effective_conduction_mass)) ** 0.5
        self.scattering_rates = cf.electron_mfp / self.speeds
        self.scattering_rates = np.column_stack((self.energies_unique, self.scattering_rates))
    
    
    def write_into_file(self):
        np.savetxt("Data/Mean travel time vs energy.csv", self.mean_travel_times, fmt='%2.4e', header="Energy [J], Travel time [s]", encoding='utf-8', delimiter=',')
        np.savetxt("Data/Transport distribution function.csv", self.transport_function, fmt='%2.4e', header="Energy [J], Transport distribution function [s^-1 J^-1 m^-1]", encoding='utf-8', delimiter=',')
        np.savetxt("Data/Electron conductivity.csv", np.column_stack((self.electron_conductivity, self.true_conductivity[:,1])), fmt='%2.4e', header="Fermi-level [J], Conductivity [S/m], Theorical Conductivity [S/m]", encoding='utf-8', delimiter=',')
        np.savetxt("Data/Seebeck coefficient.csv", np.column_stack((self.seebeck_coeff, self.true_seebeck[:,1])), fmt='%2.4e', header="Fermi-level [J], Seebeck coefficient [V/K], Theorical Seebeck coefficient [V/K]", encoding='utf-8', delimiter=',')
        np.savetxt("Data/Power factor.csv", np.column_stack((self.power_factor, self.true_power_factor[:,1])), fmt='%2.4e', header="Fermi-level [J], Power factor [W/m/K^2], Theorical Power factor [W/m/K^2]", encoding='utf-8', delimiter=',')
        np.savetxt("Data/True transport distribution function.csv", self.true_tdf, fmt='%2.4e', header="Energy [J], True TDF [s^-1 J^-1 m^-1]", encoding='utf-8', delimiter=',')
        np.savetxt("Data/Mapping constant.csv", self.mapping_constant, fmt='%2.4e', header="Fermi-level [J], Mapping constant [m^2], Constant fit [m^2]", encoding='utf-8', delimiter=',')
        np.savetxt("Data/Electron thermal conductivity.csv", np.column_stack((self.electron_thermal_conductivity, self.true_thermal_conductivity[:,1])), fmt='%2.4e', header="Fermi-level [J], Electron thermal conductivity [W/m/K], Theorical Electron thermal conductivity [W/m/K]", encoding='utf-8', delimiter=',')
        np.savetxt("Data/Figure of merit.csv", np.column_stack((self.figure_of_merit, self.true_figure_of_merit[:,1])), fmt='%2.4e', header="Fermi-level [J], Figure of merit [unitless], Theorical Figure of merit [unitless]", encoding='utf-8', delimiter=',')
        np.savetxt("Data/Scattering rate vs energy.csv", self.scattering_rates, fmt='%2.4e', header="Energy [J], Scattering rate [s]", encoding='utf-8', delimiter=',')
    
    def compute(self):
        self.compute_unique_values()
        self.compute_physical_functions()
        self.compute_mean_travel_times()
        self.compute_true_tdf()
        self.compute_mapping_constant()
        self.compute_transport_function()
        self.compute_electron_conductivity()
        self.compute_seebeck()
        self.compute_power_factor()
        self.compute_electronic_thermal_conductivity()
        self.compute_figure_of_merit()
        self.compute_scattering_rate()
