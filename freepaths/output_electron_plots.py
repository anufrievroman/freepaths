"""Plots specific to electron transport simulations."""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import electron_volt, k as k_B

from freepaths.config import cf
from freepaths.materials import get_media_class


def interpolate_property(energy_levels: np.ndarray,
                         property_values: np.ndarray,
                         fermi_level: float) -> float:
    """Interpolate the material property at a given Fermi level."""
    return np.interp(fermi_level, energy_levels, property_values)


def plot_travel_time_vs_energy():
    """Plot mean travel time vs energy on log-log with linear regression fit."""
    energy, travel_time = np.genfromtxt(
        "Data/Mean travel time vs energy.csv",
        unpack=True, delimiter=',', usecols=(0, 1), skip_header=1)

    x = energy * 1e3 / electron_volt  # Energy in meV
    y = travel_time * 1e9             # Travel time in ns

    logx = np.log(x)
    logy = np.log(y)
    m, b = np.polyfit(logx, logy, 1)
    y_fit = np.exp(b) * x**m

    fig, ax = plt.subplots()
    ax.loglog(x, y, '-o', markersize=2, c='royalblue', label='τ(E)')
    ax.loglog(x, y_fit, '-', c='green', linewidth=1.5, label=f'Linear fit: y ~ x^{m:.2f}')
    ax.set_xlabel('Energy (meV)')
    ax.set_ylabel('Mean travel time (ns)')
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.legend(loc='best')
    fig.savefig("Travel time vs energy.pdf", format="pdf", bbox_inches="tight")
    plt.close(fig)


def plot_transport_function():
    """Plot transport distribution function vs energy (Priyadarshi et al. 2023, Fig. 3)."""
    e1, tdf1 = np.genfromtxt(
        "Data/Transport distribution function.csv",
        unpack=True, delimiter=',', usecols=(0, 1), skip_header=1)
    e2, tdf2 = np.genfromtxt(
        "Data/True transport distribution function.csv",
        unpack=True, delimiter=',', usecols=(0, 1), skip_header=1)

    x  = e1 * 1e3 / electron_volt       # Energy in meV
    y1 = tdf1 * electron_volt * 1e-25   # MC TDF in 10^25 m^-1 s^-1 eV^-1
    y2 = tdf2 * electron_volt * 1e-25   # BTE TDF

    # Ξ × (-∂f/∂E): energy window contributing to conductivity (Eq. 4)
    ef = cf.media_fermi_level if cf.media_fermi_level is not None else 0.0
    eta = (e1 - ef) / (k_B * cf.temp)
    f = 1.0 / (1.0 + np.exp(np.clip(eta, -500, 500)))
    neg_df_dE = f * (1 - f) / (k_B * cf.temp)
    prod1 = tdf1 * neg_df_dE * electron_volt**2 * 1e-25   # MC  Ξ×(-∂f/∂E) in 10^25 m^-1 s^-1 eV^-2
    prod2 = tdf2 * neg_df_dE * electron_volt**2 * 1e-25   # BTE Ξ×(-∂f/∂E)

    ef_mev = ef * 1e3 / electron_volt

    import matplotlib as mpl
    _green = mpl.cm.viridis(0.55)

    fig, ax1 = plt.subplots()
    ax1.axvline(x=ef_mev, color='gray', linestyle='--', linewidth=1, zorder=1, label=f'$E_F$ = {ef_mev:.0f} meV')
    ax1.plot(x, y2, '--', c='deeppink', linewidth=0.8, label='BTE', zorder=2)
    ax1.plot(x, y1, 'o', markersize=3, c='royalblue', linewidth=0.8, label='MC', zorder=3)
    ax1.set_xlabel('Energy (meV)')
    ax1.set_ylabel('TDF, Ξ ($10^{25}$ m$^{-1}$ s$^{-1}$ eV$^{-1}$)')
    ax1.set_xlim(left=0)
    ax1.set_ylim(bottom=0)
    ax1.legend(loc='best')

    ax2 = ax1.twinx()
    ax2.plot(x, prod2, '--', c=_green, linewidth=0.8, alpha=0.5, zorder=1)
    ax2.plot(x, prod1, 'o', markersize=3, c=_green, alpha=0.5, zorder=2)
    ax2.set_ylabel('Ξ × (−$\\partial f/\\partial E$) ($10^{25}$ m$^{-1}$ s$^{-1}$ eV$^{-2}$)', color=_green)
    ax2.tick_params(axis='y', labelcolor=_green)
    ax2.set_ylim(bottom=0)

    fig.savefig("Transport distribution function.pdf", format="pdf", bbox_inches="tight")
    plt.close(fig)


def plot_electron_conductivity():
    """Plot electron conductivity vs Fermi level (Priyadarshi et al. 2023, Eq. 4)."""
    material = get_media_class(cf.media)(cf.temp, fermi_level=cf.media_fermi_level)
    fermi_levels, conductivity, true_conductivity = np.genfromtxt(
        "Data/Electron conductivity.csv", unpack=True, delimiter=',', usecols=(0, 1, 2), skip_header=1)

    material_conductivity = interpolate_property(fermi_levels, conductivity, material.fermi_level)

    fig, ax = plt.subplots()
    ax.axhline(y=material_conductivity * 1e-3, color='gray', linestyle='--', linewidth=1, zorder=1,
               label=f"σ = {material_conductivity * 1e-3:.2f} kS/m "
                     f"({material_conductivity / 100:.2f} $\\Omega^{{-1}}$cm$^{{-1}}$) "
                     f"at {material.fermi_level * 1e3 / electron_volt:.2e} meV")
    ax.axvline(x=material.fermi_level * 1e3 / electron_volt, color='gray', linestyle='--', linewidth=1, zorder=1)
    ax.plot(fermi_levels * 1e3 / electron_volt, true_conductivity * 1e-3,
            '--', c='deeppink', linewidth=0.8, label="BTE", zorder=2)
    ax.plot(fermi_levels * 1e3 / electron_volt, conductivity * 1e-3,
            'o', markersize=3, c='royalblue', linewidth=0.8, label="MC", zorder=3)
    ax.set_xlabel('Fermi level (meV)')
    ax.set_ylabel('Electron conductivity (kS/m)')
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.legend()
    fig.savefig("Electron conductivity.pdf", format="pdf", bbox_inches="tight")
    plt.close(fig)


def plot_seebeck_coefficient():
    """Plot Seebeck coefficient vs Fermi level (Priyadarshi et al. 2023, Eq. 5)."""
    material = get_media_class(cf.media)(cf.temp, fermi_level=cf.media_fermi_level)
    fermi_levels, seebeck, true_seebeck = np.genfromtxt(
        "Data/Seebeck coefficient.csv", unpack=True, delimiter=',', usecols=(0, 1, 2), skip_header=1)

    material_seebeck = abs(interpolate_property(fermi_levels, seebeck, material.fermi_level) * 1e3)

    fig, ax = plt.subplots()
    ax.axhline(y=material_seebeck, color='gray', linestyle='--', linewidth=1, zorder=1,
               label=f"|S| = {material_seebeck:.2f} mV/K at {material.fermi_level * 1e3 / electron_volt:.2e} meV")
    ax.axvline(x=material.fermi_level * 1e3 / electron_volt, color='gray', linestyle='--', linewidth=1, zorder=1)
    ax.plot(fermi_levels * 1e3 / electron_volt, np.abs(true_seebeck) * 1e3,
            '--', c='deeppink', linewidth=0.8, label="BTE", zorder=2)
    ax.plot(fermi_levels * 1e3 / electron_volt, np.abs(seebeck) * 1e3,
            'o', markersize=3, c='royalblue', linewidth=0.8, label="MC", zorder=3)
    ax.set_xlabel('Fermi level (meV)')
    ax.set_ylabel('|Seebeck coefficient| (mV/K)')
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.legend()
    fig.savefig("Seebeck coefficient.pdf", format="pdf", bbox_inches="tight")
    plt.close(fig)


def plot_power_factor():
    """Plot thermoelectric power factor vs Fermi level (Priyadarshi et al. 2023, Eq. 6)."""
    material = get_media_class(cf.media)(cf.temp, fermi_level=cf.media_fermi_level)
    fermi_levels, power_factor, true_power_factor = np.genfromtxt(
        "Data/Power factor.csv", unpack=True, delimiter=',', usecols=(0, 1, 2), skip_header=1)

    material_pf = interpolate_property(fermi_levels, power_factor, material.fermi_level) * 1e3

    fig, ax = plt.subplots()
    ax.axhline(y=material_pf, color='gray', linestyle='--', linewidth=1, zorder=1,
               label=f"PF = {material_pf:.2f} mW/m·K² at {material.fermi_level * 1e3 / electron_volt:.2e} meV")
    ax.axvline(x=material.fermi_level * 1e3 / electron_volt, color='gray', linestyle='--', linewidth=1, zorder=1)
    ax.plot(fermi_levels * 1e3 / electron_volt, true_power_factor * 1e3,
            '--', c='deeppink', linewidth=0.8, label="BTE", zorder=2)
    ax.plot(fermi_levels * 1e3 / electron_volt, power_factor * 1e3,
            'o', markersize=3, c='royalblue', linewidth=0.8, label="MC", zorder=3)
    ax.set_xlabel('Fermi level (meV)')
    ax.set_ylabel('Power factor (mW/m·$K^{2}$)')
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.legend()
    fig.savefig("Power factor.pdf", format="pdf", bbox_inches="tight")
    plt.close(fig)


def plot_mapping_constant():
    """Plot the mapping constant C vs Fermi level."""
    material = get_media_class(cf.media)(cf.temp, fermi_level=cf.media_fermi_level)
    fermi_level, mapping_constant = np.genfromtxt(
        "Data/Mapping constant.csv", unpack=True, delimiter=',', usecols=(0, 1), skip_header=1)

    mean_constant = np.mean(mapping_constant)

    fig, ax = plt.subplots()
    ax.axhline(y=mean_constant, color='gray', linestyle='--', linewidth=1,
               label=f"$\\langle C \\rangle$ = {mean_constant:.4e} m²")
    ax.plot(fermi_level * 1e3 / electron_volt, mapping_constant, '-o', markersize=2, c='royalblue', label='MC')
    ax.set_xlabel('Fermi level (meV)')
    ax.set_ylabel('Mapping constant (m²)')
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.legend()
    fig.savefig("Mapping constant.pdf", format="pdf", bbox_inches="tight")
    plt.close(fig)


def plot_electron_thermal_conductivity():
    """Plot electronic thermal conductivity vs Fermi level (Priyadarshi et al. 2023, Eq. 7)."""
    material = get_media_class(cf.media)(cf.temp, fermi_level=cf.media_fermi_level)
    fermi_level, thermal_conductivity, true_thermal_conductivity = np.genfromtxt(
        "Data/Electron thermal conductivity.csv", unpack=True, delimiter=',', usecols=(0, 1, 2), skip_header=1)

    material_kappa = interpolate_property(fermi_level, thermal_conductivity, material.fermi_level)

    fig, ax = plt.subplots()
    ax.axhline(y=material_kappa, color='gray', linestyle='--', linewidth=1,
               label=f"κ_e = {material_kappa:.4f} W/m·K at {material.fermi_level * 1e3 / electron_volt:.2e} meV")
    ax.axvline(x=material.fermi_level * 1e3 / electron_volt, color='gray', linestyle='--', linewidth=1)
    ax.plot(fermi_level * 1e3 / electron_volt, true_thermal_conductivity,
            '--', c='deeppink', linewidth=0.8, label="BTE", zorder=2)
    ax.plot(fermi_level * 1e3 / electron_volt, thermal_conductivity,
            'o', markersize=3, c='royalblue', linewidth=0.8, label="MC", zorder=3)
    ax.set_xlabel('Fermi level (meV)')
    ax.set_ylabel('Thermal conductivity (W/m·K)')
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.legend()
    fig.savefig("Electron thermal conductivity.pdf", format="pdf", bbox_inches="tight")
    plt.close(fig)


def plot_scattering_rate_vs_energy():
    """Plot internal scattering relaxation time τ(E) = mfp/v(E) vs energy."""
    energies, scattering_rate = np.genfromtxt(
        "Data/Scattering rate vs energy.csv", unpack=True, delimiter=',', usecols=(0, 1), skip_header=1)

    fig, ax = plt.subplots()
    ax.plot(energies * 1e3 / electron_volt, scattering_rate * 1e12, '-o', markersize=2, c='royalblue')
    ax.set_xlabel('Energy (meV)')
    ax.set_ylabel('Relaxation time (ps)')
    fig.savefig("Scattering rate vs energy.pdf", format="pdf", bbox_inches="tight")
    plt.close(fig)
