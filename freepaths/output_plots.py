"""Module that calculates and outputs vaious plots and distributions from the saved files"""

import logging
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LogNorm
from matplotlib import font_manager
from scipy.constants import electron_volt

from freepaths.config import cf
from freepaths.particle_types import ParticleType
from freepaths.materials import get_media_class
from freepaths.output_structure import draw_structure_top_view, draw_structure_side_view
from freepaths.materials import Si, SiC, Graphite, Ge
import matplotlib.pyplot as plt
import matplotlib as mpl 
mpl.rcParams['pdf.compression'] = 9   # compresse PDF flux


# Style of the plots:
all_fonts = font_manager.get_font_names()
if 'Arial' in all_fonts:
    plt.rcParams['font.family'] = 'Arial'
else:
    # Use a generic font or specify multiple options
    logging.warning("Arial font not available. Falling back on default sans-serif font")
    plt.rcParams['font.family'] = ['sans-serif']
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.titlesize'] = 8
plt.rcParams['lines.linewidth'] = 1.0
plt.rcParams['xtick.direction'] = "in"
plt.rcParams['ytick.direction'] = "in"
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['xtick.minor.visible'] = False
plt.rcParams['ytick.minor.visible'] = False
plt.rcParams['legend.frameon'] = False
plt.rcParams['figure.autolayout'] = True
plt.rcParams['figure.figsize'] = [5, 3.5]
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['grid.linewidth'] = 0.5


def distribution_calculation(filename, data_range, number_of_nodes):
    """Calculate distribution of numbers (histogram) in a given file"""
    data = np.loadtxt(filename, encoding='utf-8')
    if data_range is None:
        data_range = np.max(data)
    distribution = np.zeros((number_of_nodes, 2))
    distribution[:, 0] = np.linspace(0, data_range, number_of_nodes)
    distribution[:, 1], _ = np.histogram(data[data != 0], number_of_nodes, range=(0, data_range))
    return distribution

# --- Helper: binning 2D avec moyenne de la transmission ---

def _bin_average_2d(x, y, z, nx=200, ny=100, min_count=1, x_range=None, y_range=None): #add 22/08
    """
    Moyenne z dans une grille (x,y).
    Retourne des points agrégés (x_c, y_c), z_mean et le count par bin.
    """
    x = np.asarray(x); y = np.asarray(y); z = np.asarray(z)

    if x_range is None: x_range = (np.nanmin(x), np.nanmax(x))
    if y_range is None: y_range = (np.nanmin(y), np.nanmax(y))

    xb = np.linspace(x_range[0], x_range[1], nx + 1)
    yb = np.linspace(y_range[0], y_range[1], ny + 1)

    ix = np.digitize(x, xb) - 1
    iy = np.digitize(y, yb) - 1

    mask = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny) & np.isfinite(z)
    ix = ix[mask]; iy = iy[mask]; z = z[mask]

    sumz = np.zeros((nx, ny), dtype=float)
    cnt  = np.zeros((nx, ny), dtype=int)
    np.add.at(sumz, (ix, iy), z)
    np.add.at(cnt,  (ix, iy), 1)

    mean = np.zeros_like(sumz); nonzero = cnt > 0
    mean[nonzero] = sumz[nonzero] / cnt[nonzero]

    xc = 0.5 * (xb[:-1] + xb[1:])
    yc = 0.5 * (yb[:-1] + yb[1:])

    I, J = np.where(cnt >= min_count)
    return xc[I], yc[J], mean[I, J], cnt[I, J]

# --- Helper: binning 1D avec moyenne ---
def _bin_average_1d(x, y, nbins=120, min_count=1, x_range=None): # add 22/08
    """
    Moyenne y dans des bins de x.
    Retourne: x_centers, mean(y), std(y), count par bin (filtré par min_count).
    """
    x = np.asarray(x); y = np.asarray(y)

    if x_range is None:
        x_range = (np.nanmin(x), np.nanmax(x))

    edges = np.linspace(x_range[0], x_range[1], nbins + 1)
    idx = np.digitize(x, edges) - 1

    mask = (idx >= 0) & (idx < nbins) & np.isfinite(y)
    idx = idx[mask]; y = y[mask]

    sumy  = np.zeros(nbins); sumy2 = np.zeros(nbins); cnt = np.zeros(nbins, dtype=int)
    np.add.at(sumy,  idx, y)
    np.add.at(sumy2, idx, y*y)
    np.add.at(cnt,   idx, 1)

    mean = np.full(nbins, np.nan); std = np.full(nbins, np.nan)
    nz = cnt > 0
    mean[nz] = sumy[nz] / cnt[nz]
    var = np.maximum(sumy2[nz] / cnt[nz] - mean[nz]**2, 0.0)
    std[nz]  = np.sqrt(var)

    centers = 0.5 * (edges[:-1] + edges[1:])
    keep = cnt >= min_count
    return centers[keep], mean[keep], std[keep], cnt[keep]



def angle_distribution_calculation():
    """Analyse measured particle angles and create their distribution"""
    all_exit_angles = np.loadtxt("Data/All exit angles.csv", dtype='float', encoding='utf-8')
    initial_angles = np.loadtxt("Data/All initial angles.csv", dtype='float', encoding='utf-8')
    hole_diff_angles = np.loadtxt("Data/All hole diffuse scattering angles.csv", dtype='float', encoding='utf-8')
    hole_spec_angles = np.loadtxt("Data/All hole specular scattering angles.csv", dtype='float', encoding='utf-8')
    distribution = np.zeros((360, 5)) 
    distribution[:, 0] = range(-180, 180)
    exit_angles = all_exit_angles[all_exit_angles != 0]
    distribution[:, 1], _ = np.histogram(np.degrees(exit_angles), 360, range=(-180, 180))
    distribution[:, 2], _ = np.histogram(np.degrees(initial_angles), 360, range=(-180, 180))
    distribution[:, 3], _ = np.histogram(np.degrees(hole_diff_angles), 360, range=(-180, 180)) 
    distribution[:, 4], _ = np.histogram(np.degrees(hole_spec_angles), 360, range=(-180, 180)) 
    return distribution

def interfaces_transmission_angles_calculation(): 
    all_exit_angles = np.loadtxt("Data/All exit angles.csv", dtype='float', encoding='utf-8')
    initial_angles = np.loadtxt("Data/All initial angles.csv", dtype='float', encoding='utf-8')
    interfaces_transmission_specular_angles = np.loadtxt("Data/All interfaces transmission specular.csv", dtype='float', encoding='utf-8')  
    interfaces_transmission_diffuse_angles = np.loadtxt("Data/All interfaces transmission diffuse.csv", dtype='float', encoding='utf-8')  
    interfaces_angles = np.loadtxt("Data/All interfaces angles.csv", dtype='float', encoding='utf-8')  
    distribution = np.zeros((360, 6)) 
    distribution[:, 0] = range(-180, 180)
    exit_angles = all_exit_angles[all_exit_angles != 0]
    distribution[:, 1], _ = np.histogram(np.degrees(exit_angles), 360, range=(-180, 180))
    distribution[:, 2], _ = np.histogram(np.degrees(initial_angles), 360, range=(-180, 180))
    distribution[:, 3], _ = np.histogram(np.degrees(interfaces_transmission_specular_angles), 360, range=(-180, 180))  
    distribution[:, 4], _ = np.histogram(np.degrees(interfaces_transmission_diffuse_angles), 360, range=(-180, 180))  
    distribution[:, 5], _ = np.histogram(np.degrees(interfaces_angles), 360, range=(-180, 180))  
    return distribution


def scattering_angle_distribution_calculation():
    """Analyse scattering particle angles and create their distribution"""
    hole_diff_angles = np.loadtxt("Data/All hole diffuse scattering angles.csv", dtype='float', encoding='utf-8')
    hole_spec_angles = np.loadtxt("Data/All hole specular scattering angles.csv", dtype='float', encoding='utf-8')
    distribution = np.zeros((360, 3)) 
    distribution[:, 0] = range(-180, 180) 
    distribution[:, 1], _ = np.histogram(np.degrees(hole_diff_angles), 360, range=(-180, 180))
    distribution[:, 2], _ = np.histogram(np.degrees(hole_spec_angles), 360, range=(-180, 180))
    return distribution

def scattering_interfaces_angles_distribution_calculation(): 
    interfaces_transmission_specular_angles = np.loadtxt("Data/All interfaces transmission specular.csv", dtype='float', encoding='utf-8') 
    interfaces_transmission_diffuse_angles = np.loadtxt("Data/All interfaces transmission diffuse.csv", dtype='float', encoding='utf-8') 
    interfaces_angles = np.loadtxt("Data/All interfaces angles.csv", dtype='float', encoding='utf-8') 
    distribution = np.zeros((360, 4)) 
    distribution[:, 0] = range(-180, 180) 
    distribution[:, 1], _ = np.histogram(np.degrees(interfaces_transmission_specular_angles), 360, range=(-180, 180)) 
    distribution[:, 2], _ = np.histogram(np.degrees(interfaces_transmission_diffuse_angles), 360, range=(-180, 180)) 
    distribution[:, 3], _ = np.histogram(np.degrees(interfaces_angles), 360, range=(-180, 180)) 
    return distribution

def wavelength_distribution_calculation(number_of_nodes):
    """Calculate particle wavelength distribution from their frequencies and velocities"""
    frequencies = np.loadtxt("Data/All initial frequencies.csv", encoding='utf-8')
    speeds = np.loadtxt("Data/All group velocities.csv", encoding='utf-8')
    wavelengths = np.zeros((len(speeds)))
    wavelengths[:] = speeds[:] / frequencies[:]
    data_range = np.amax(wavelengths)
    distribution = np.zeros((number_of_nodes, 2))
    distribution[:, 0] = np.linspace(0, data_range, number_of_nodes)
    distribution[:, 1], _ = np.histogram(wavelengths, number_of_nodes, range=(0, data_range))
    return distribution


def cumulative_conductivity_calculation():
    """Calculate cumulative function of thermal conductivity vs mean freee path"""
    mfp = np.loadtxt("Data/All mean free paths.csv", encoding='utf-8')
    kappa = np.loadtxt("Data/All thermal conductivities.csv", encoding='utf-8')
    sorted_indices = np.argsort(mfp)
    sorted_mfp = mfp[sorted_indices]
    sorted_kappa = kappa[sorted_indices]
    cumulative_thermal_conductivity = np.cumsum(sorted_kappa)
    return sorted_mfp, cumulative_thermal_conductivity

def interpolate_property(energy_levels: np.ndarray,
                         property_values: np.ndarray,
                         fermi_level: float) -> float:
    """
    Interpolate the material property at a given Fermi level.
    """
    return np.interp(fermi_level, energy_levels, property_values)

def plot_cumulative_thermal_conductivity(mfp_sampling):
    """Plot distribution cumulative thermal conductivity vs mean free path"""
    if not mfp_sampling:
        return
    mfp, kappa = cumulative_conductivity_calculation()
    fig, ax = plt.subplots()
    ax.plot(mfp * 1e6, kappa, 'royalblue')
    ax.set_xlabel('Mean free path (μm)')
    ax.set_ylabel('Cumulative thermal conductivity (W/m·K)')
    ax.set_xscale('log')
    fig.savefig("Distribution of thermal conductivity.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)
    np.savetxt('Data/Distribution of thermal conductivity.csv', np.vstack((mfp, kappa)).T, fmt='%1.3e', delimiter=",")


def plot_angle_distribution():
    """Plot distribution of initial and exit angles"""
    angle_distributions = angle_distribution_calculation()
    fig, ax = plt.subplots()
    ax.plot(angle_distributions[:, 0], angle_distributions[:, 1], 'royalblue')
    ax.plot(angle_distributions[:, 0], angle_distributions[:, 2], 'deeppink')
    ax.set_xlabel('Angle (degree)')
    ax.set_ylabel('Number of particles')
    ax.set_ylim(bottom=0)
    ax.legend(["At cold side", "At hot side"])
    fig.savefig("Distribution of angles.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)
    np.savetxt('Data/Distribution of angles.csv', angle_distributions, fmt='%1.3e', delimiter=",")


def plot_scattering_angle_distribution():
    """Plot distribution of hole scattering angles"""
    angle_distributions = scattering_angle_distribution_calculation()
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.grid(zorder=0)
    ax.plot(np.deg2rad(angle_distributions[:, 0]), angle_distributions[:, 1], 'royalblue', label="Diffuse", zorder=2)
    ax.fill_between(np.deg2rad(angle_distributions[:, 0]), angle_distributions[:, 1], 0, alpha=0.2, zorder=2)
    ax.plot(np.deg2rad(angle_distributions[:, 0]), angle_distributions[:, 2], 'deeppink', label="Specular", zorder=3)
    ax.fill_between(np.deg2rad(angle_distributions[:, 0]), angle_distributions[:, 2], 0, alpha=0.2, zorder=3)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.legend(facecolor='white', framealpha=1, ncols=2, loc="lower center", bbox_to_anchor=(0.5, -0.17))
    ax.set_title('Number of scattered particles per angle')
    fig.savefig("Distribution of hole scattering angles.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)
    np.savetxt('Data/Distribution of hole scattering angles.csv', angle_distributions, fmt='%1.3e', delimiter=",")

def plot_interfaces_angles_distribution(): 
    """Plot distribution of interfaces angles"""   
    angle_distributions = scattering_interfaces_angles_distribution_calculation()
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.grid(zorder=0)
    ax.plot(np.deg2rad(angle_distributions[:, 0]), angle_distributions[:, 1], 'royalblue', label="Diffuse", zorder=2)  
    ax.fill_between(np.deg2rad(angle_distributions[:, 0]), angle_distributions[:, 1], 0, alpha=0.2, zorder=2)  
    ax.plot(np.deg2rad(angle_distributions[:, 0]), angle_distributions[:, 2], 'deeppink', label="Specular", zorder=3)  
    ax.fill_between(np.deg2rad(angle_distributions[:, 0]), angle_distributions[:, 2], 0, alpha=0.2, zorder=3)  
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.legend(facecolor='white', framealpha=1, ncols=2, loc="lower center", bbox_to_anchor=(0.5, -0.17))
    ax.set_title('Number of transmitted phonons per angle')
    fig.savefig("Distribution of interfaces angles.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)
    np.savetxt('Data/Distribution of interfaces angles.csv', angle_distributions, fmt='%1.3e', delimiter=",")


import subprocess
import os

def open_pdf_with_mupdf(pdf_path):
    try:
        # path of mu pdf
        mupdf_path = r"C:\Users\victo\Downloads\mupdf-1.26.2-windows\mupdf-1.26.2-windows\mupdf-gl.exe"

        # Lancer mupdf avec ton fichier
        subprocess.Popen([mupdf_path, os.path.abspath(pdf_path)])
    except Exception as e:
        print(f"[ERROR] Impossible d'ouvrir le PDF avec MuPDF : {e}")

def plot_transmission_vs_angle():
    try:
        angles = np.loadtxt("Data/All interfaces angles.csv", dtype='float', encoding='utf-8')
        angles_degrees = np.degrees(angles)
        transmission_factor = np.loadtxt("Data/All interfaces transmission factor.csv", dtype='float', encoding='utf-8')
        modes = np.loadtxt("Data/All interfaces mode.csv", dtype='int', encoding='utf-8')

        if not (len(angles) == len(transmission_factor) == len(modes)):
            raise ValueError("Data arrays (angles, transmission, modes) must have the same length.")

        fig, ax = plt.subplots(figsize=(8, 5))
        mode_colors = {0: 'red', 1: 'blue', 2: 'green'}
        mode_labels = {0: 'LA', 1: 'TA1', 2: 'TA2'}

        for mode_number in [0, 1, 2]:
            mask = (modes == mode_number)
            if np.any(mask):
                ax.scatter(
                    angles_degrees[mask],
                    transmission_factor[mask],
                    alpha=0.5,
                    s=10,
                    c=mode_colors[mode_number],
                    label=mode_labels[mode_number],
                    edgecolors='none',
                    rasterized=True,          
                )

        ax.set_xlabel('Angle (degree)')
        ax.set_ylabel('Transmission factor')
        ax.set_title('Transmission factor vs Angle (colored by mode)')
        ax.legend()
        ax.set_xlim(0, 90)
        ax.set_ylim(0, 1)
        ax.grid(True, linestyle='--', alpha=0.5)
        fig.tight_layout()

        out_pdf = "Transmission coefficient vs angle.pdf"
        fig.savefig(out_pdf, dpi=220, format='pdf', bbox_inches="tight")  # dpi = résolution
        plt.close(fig)

    except Exception as e:
        print(f"[ERROR] plot_transmission_vs_angle: {e}")

def plot_transmission_vs_wavelength():
    """Plot transmission vs wavelength"""
    wavelength = np.loadtxt("Data/All interfaces wavelength.csv", dtype='float', encoding='utf-8')
    transmission_factor = np.loadtxt("Data/All interfaces transmission factor.csv", dtype='float', encoding='utf-8')

    # check the size correspondance and if there are data 
    if len(wavelength) == 0 or len(wavelength) != len(transmission_factor):
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    
    # convert wavelength to nm for better readability
    ax.scatter(wavelength * 1e9, transmission_factor, alpha=0.7, s=10)

    # labels
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Transmission factor')

    # axes limits (nm)
    wl_min_nm = np.min(wavelength) * 1e9
    wl_max_nm = np.max(wavelength) * 1e9
    if wl_min_nm == wl_max_nm:
        # if all lambda are identical
        if wl_min_nm == 0:
            ax.set_xlim(0, 100) # random plage
        else:
            margin = wl_min_nm * 0.1
            ax.set_xlim(wl_min_nm - margin, wl_max_nm + margin)
    else:
         ax.set_xlim(0.8 * wl_min_nm, 1.2 * wl_max_nm)

    ax.set_ylim(0, 1)
    ax.grid(True, linestyle='--', alpha=0.5)
    fig.tight_layout()
    fig.savefig("Transmission factor vs wavelength.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)

    # save the brut data 
    output_data = np.column_stack((wavelength * 1e9, transmission_factor))
    header = "Wavelength_nm,Transmission_Factor"
    np.savetxt('Data/Distribution of interfaces wavelength.csv', output_data, fmt='%.6f,%.6f', delimiter=",", header=header, comments='')

def plot_transmission_vs_frequency():
    """Plot transmission vs frequency"""
    frequency = np.loadtxt("Data/All interfaces frequency.csv", dtype='float', encoding='utf-8')
    transmission_factor = np.loadtxt("Data/All interfaces transmission factor.csv", dtype='float', encoding='utf-8')

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(frequency, transmission_factor, alpha=0.5, s=1, c='blue', label='Événements')
    ax.set_xlabel('Frequency (THz)')
    ax.set_ylabel('Transmission factor')
    ax.legend()
    # axes limits (en nm)
    wl_min = np.min(frequency) 
    wl_max = np.max(frequency) 
    if wl_min == wl_max:
        # if all lambda are identical
        if wl_min == 0:
            ax.set_xlim(0, 100) # random plage
        else:
            margin = wl_min * 0.1
            ax.set_xlim(wl_min - margin, wl_max + margin)
    else:
         ax.set_xlim(0.8 * wl_min, 1.2 * wl_max)
   
    ax.set_ylim(0, 1)
    ax.grid(True, linestyle='--', alpha=0.5)
    fig.tight_layout()
    fig.savefig("Transmission coefficient vs frequency.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)
    np.savetxt('Data/Distribution of frequency.csv', frequency, fmt='%1.3e', delimiter=",")

def plot_transmission_heatmaps_by_mode():
    """
    Génère un PDF multipage : 1 page par mode (LA, TA1, TA2).
    Les nuages sont rasterisés pour un PDF léger.
    Ouvre automatiquement le PDF avec MuPDF à la fin.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    # --- loading data ---
    freq  = np.loadtxt("Data/All interfaces frequency.csv", dtype=float, encoding="utf-8")
    ang   = np.loadtxt("Data/All interfaces angles.csv",    dtype=float, encoding="utf-8")
    T     = np.loadtxt("Data/All interfaces transmission factor.csv", dtype=float, encoding="utf-8")
    modes = np.loadtxt("Data/All interfaces mode.csv",      dtype=int,   encoding="utf-8")

    if not (len(freq) == len(ang) == len(T) == len(modes)):
        raise ValueError("The files does not have the same length.")

    ang_deg = np.degrees(ang)

    x_min, x_max = np.nanmin(freq), np.nanmax(freq)
    y_min, y_max = 0.0, 90.0
    vmin, vmax = 0.0, 1.0  # transmission e [0,1]

    mode_labels = {0: "LA", 1: "TA1", 2: "TA2"}
    pdf_filename = "Transmission_heatmaps_all_modes.pdf"

    with PdfPages(pdf_filename) as pdf:
        for m in (0, 1, 2):
            mask = (modes == m) & np.isfinite(freq) & np.isfinite(ang_deg) & np.isfinite(T)

            fig, ax = plt.subplots(figsize=(8, 6))
            if np.any(mask):
                sc = ax.scatter(
                    freq[mask], ang_deg[mask],
                    c=T[mask], cmap="jet", vmin=vmin, vmax=vmax,
                    s=10, edgecolors="none", rasterized=True  # << RASTER
                )
                cbar = plt.colorbar(sc, ax=ax)
                cbar.set_label("Transmission factor")
            else:
                ax.text(0.5, 0.5, "No data for this mode", ha="center", va="center", transform=ax.transAxes)

            ax.set_xlabel("Frequency (THz)")
            ax.set_ylabel("Incident angle (degrees)")
            ax.set_title(f"Transmission vs Frequency and Angle — Mode {mode_labels[m]}")
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(y_min, y_max)
            fig.tight_layout()

            # one page per mode
            pdf.savefig(fig, dpi=220, bbox_inches="tight")
            plt.close(fig)


def plot_free_path_distribution():
    """Plot distribution of free path"""
    filename = "Data/All free paths.csv"
    free_path_distribution = distribution_calculation(filename, None, cf.number_of_nodes)
    fig, ax = plt.subplots()
    ax.plot(free_path_distribution[:, 0] * 1e6, free_path_distribution[:, 1], 'royalblue')
    ax.set_xscale('log')
    ax.set_xlabel('Free flights (μm)')
    ax.set_ylabel('Number of flights')
    fig.savefig("Distribution of free paths.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)
    np.savetxt('Data/Distribution of free paths.csv', free_path_distribution, fmt='%1.3e', delimiter=",")


def plot_frequency_distribution():
    """Plot distribution of frequencies"""
    filename = "Data/All initial frequencies.csv"
    frequency_distribution = distribution_calculation(filename, None, cf.number_of_nodes)
    fig, ax = plt.subplots()
    ax.plot(frequency_distribution[:, 0] * 1e-12, frequency_distribution[:, 1], 'royalblue')
    ax.set_xlabel('Frequency (THz)')
    ax.set_ylabel('Number of particles')
    fig.savefig("Distribution of initial frequencies.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)
    np.savetxt('Data/Distribution of initial frequencies.csv', frequency_distribution, fmt='%1.3e', delimiter=",")


def plot_wavelength_distribution():
    """Plot distribution of wavelength"""
    wavelength_distribution = wavelength_distribution_calculation(cf.number_of_nodes)
    fig, ax = plt.subplots()
    ax.plot(wavelength_distribution[:, 0] * 1e9, wavelength_distribution[:, 1], 'royalblue')
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Number of particles')
    fig.savefig("Distribution of wavelengths.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)
    np.savetxt('Data/Distribution of wavelengths.csv', wavelength_distribution, fmt='%1.3e', delimiter=",")


def plot_travel_time_distribution():
    """Plot distribution of wavelength"""
    travel_time_distribution = distribution_calculation("Data/All travel times.csv", None, cf.number_of_nodes)
    fig, ax = plt.subplots()
    ax.plot(travel_time_distribution[:, 0] * 1e9, travel_time_distribution[:, 1], 'royalblue')
    ax.set_xscale('log')
    ax.set_xlabel('Travel time (ns)')
    ax.set_ylabel('Number of particles')
    fig.savefig("Distribution of travel times.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)
    np.savetxt('Data/Distribution of travel times.csv', travel_time_distribution, fmt='%1.3e', delimiter=",")


def plot_mean_free_path_distribution():
    """Plot distribution of MFP per particle"""
    mean_free_path_distribution = distribution_calculation("Data/All mean free paths.csv", None, cf.number_of_nodes)
    fig, ax = plt.subplots()
    ax.plot(mean_free_path_distribution[:, 0] * 1e9, mean_free_path_distribution[:, 1], 'royalblue')
    ax.set_xscale('log')
    ax.set_xlabel('Mean free path (nm)')
    ax.set_ylabel('Number of particles')
    fig.savefig("Distribution of MFPs.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)
    np.savetxt('Data/Distribution of MFPs.csv', mean_free_path_distribution, fmt='%1.3e', delimiter=",")


def plot_velocity_distribution():
    """Plot distribution of group velocities"""
    fig, ax = plt.subplots()
    speeds = np.loadtxt("Data/All group velocities.csv")
    frequencies = np.loadtxt("Data/All initial frequencies.csv")
    ax.plot(frequencies, speeds, '.', markersize=2, c='royalblue')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Group velocity (m/s)')
    fig.savefig('Group velocities.pdf', format='pdf', bbox_inches="tight")
    plt.close(fig)

def plot_energy_distribution():
    """Plot distribution of energy"""
    fig, ax = plt.subplots()
    energy_distribution = distribution_calculation("Data/All initial energies.csv", None, cf.number_of_nodes)
    ax.plot(energy_distribution[:, 0] / (electron_volt*1e-3), energy_distribution[:, 1], 'royalblue') # Convert J to meV
    ax.set_xlabel('Energy (meV)')
    ax.set_ylabel('Number of particles')
    fig.savefig("Distribution of energy.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)

def plot_travel_time_vs_energy():
    """Plot mean travel time vs energy, slope on log-log and linear regression fit."""
    # Load data
    energy, travel_time = np.genfromtxt(
        "Data/Mean travel time vs energy.csv",
        unpack=True,
        delimiter=',',
        usecols=(0, 1),
        skip_header=1
    )

    # Convert to plotting units
    x = energy * 1e3 / electron_volt  # Energy in meV
    y = travel_time * 1e9             # Travel time in ns

    # Compute log values
    logx = np.log(x)
    logy = np.log(y)

    # Linear regression in log-log space: fit log(y) = m*log(x) + b
    m, b = np.polyfit(logx, logy, 1)
    # Regression line in original scale: y_fit = exp(b) * x**m
    y_fit = np.exp(b) * x**m


    # Set up figure and primary axis
    fig, ax1 = plt.subplots()
    ax1.loglog(x, y, '-o', markersize=2, c='royalblue', label='τ(E)')
    ax1.loglog(x, y_fit, '-', c='green', linewidth=1.5,
               label=f'Regression: y as x^{m:.2f}')
    ax1.set_xlabel('Energy (meV)')
    ax1.set_ylabel('Mean travel time (ns)', color='royalblue')
    ax1.grid(True, linestyle='--', alpha=0.7)

    # Combine legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    ax1.legend(lines1, labels1, loc='best')

    # Save and close
    fig.savefig("Travel time vs energy.pdf", format="pdf", bbox_inches="tight")
    plt.close(fig)

def plot_transport_function():
    """Plot transport distribution function vs energy with regression and slope."""
    # Load computed and theoretical TDF
    e1, tdf1 = np.genfromtxt(
        "Data/Transport distribution function.csv",
        unpack=True, delimiter=',', usecols=(0,1), skip_header=1)
    e2, tdf2 = np.genfromtxt(
        "Data/True transport distribution function.csv",
        unpack=True, delimiter=',', usecols=(0,1), skip_header=1)

    # Convert units
    x = e1 * 1e3 / electron_volt                   # Energy in meV
    y1 = tdf1 * electron_volt * 1e-25             # Computed TDF in 10^25 m^-1 s^-1 eV^-1
    y2 = tdf2 * electron_volt * 1e-25             # Theoretical TDF

    # Log variables for computed TDF
    logx = np.log(x)
    logy1 = np.log(y1)

    # Regression log(y1) = m*log(x) + b
    m, b = np.polyfit(logx, logy1, 1)
    y1_fit = np.exp(b) * x**m


    # Plotting
    fig, ax1 = plt.subplots()
    ax1.loglog(x, y1, '-o', markersize=2, c='royalblue', label='Computed TDF')
    ax1.loglog(x, y1_fit, '-', c='green', linewidth=1.5,
               label=f'Regression: TDF as E^{m:.2f}')
    ax1.loglog(x, y2, '-o', markersize=2, c='darkorange', label='Theoretical TDF')
    ax1.set_xlabel('Energy (meV)')
    ax1.set_ylabel('Transport distribution function ($10^{25}$ m$^{-1}$ s$^{-1}$ eV$^{-1}$)', color='black')
    ax1.grid(True, linestyle='--', alpha=0.7)


    # Legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    ax1.legend(lines1, labels1, loc='best')

    # Save
    fig.savefig("Transport distribution function.pdf", format="pdf", bbox_inches="tight")
    plt.close(fig)

def plot_electron_conductivity():
    """Plot electron conductivity with respect to fermi-level"""
    material = get_media_class(cf.media)(cf.temp, fermi_level=cf.media_fermi_level)
    fig, ax = plt.subplots()
    fermi_levels, conductivity, theorical_conductivity = np.genfromtxt("Data/Electron conductivity.csv", unpack=True, delimiter=',', usecols=(0,1,2), skip_header=1)
    ax.plot(fermi_levels * 1e3 / electron_volt, conductivity, '-o', markersize=2, c='royalblue', label="Computed")
    # ax.plot(fermi_levels * 1e3 / electron_volt, theorical_conductivity, '-o', markersize=2, c='darkorange', label="Bulk/Analytical")
    ax.set_xlabel('Fermi-level (meV)')
    ax.set_ylabel('Electron conductivity (S/m)')
    ax.grid(True, linestyle='--', alpha=0.7)
    material_conductivity = interpolate_property(fermi_levels, conductivity, material.fermi_level)
    ax.axhline(
        y=material_conductivity,
        color='gray',
        linestyle='--',
        linewidth=1,
        label=f"y = {material_conductivity:.4e}"
    )
    ax.axvline(
        x=material.fermi_level * 1e3 / electron_volt,
        color='gray',
        linestyle='--',
        linewidth=1,
        label=f"x = {material.fermi_level*1e3 / electron_volt:.2e}"
    )
    plt.legend()
    fig.savefig("Electron conductivity.pdf", format="pdf", bbox_inches="tight")
    plt.close(fig)

def plot_seebeck_coefficient():
    """Plot the Seebeck coefficient with respect to fermi-level"""
    material = get_media_class(cf.media)(cf.temp, fermi_level=cf.media_fermi_level)
    fig, ax = plt.subplots()
    fermi_levels, seebeck, theorical_seebeck = np.genfromtxt("Data/Seebeck coefficient.csv", unpack=True, delimiter=',', usecols=(0,1,2), skip_header=1)
    ax.plot(fermi_levels * 1e3 / electron_volt, seebeck * 1e3, '-o', markersize=2, c='royalblue', label="Computed")
    # ax.plot(fermi_levels * 1e3 / electron_volt, theorical_seebeck * 1e3, '-o', markersize=2, c='darkorange', label="Bulk/Analytical")
    ax.set_xlabel('Fermi-level (meV)')
    ax.set_ylabel('Seebeck coefficient (mV/K)')
    ax.grid(True, linestyle='--', alpha=0.7)
    material_seebeck = interpolate_property(fermi_levels, seebeck, material.fermi_level) * 1e3
    ax.axhline(
        y=material_seebeck,
        color='gray',
        linestyle='--',
        linewidth=1,
        label=f"y = {material_seebeck:.2e}"
    )
    plt.legend()
    fig.savefig("Seebeck coefficient.pdf", format="pdf", bbox_inches="tight")
    plt.close(fig)

def plot_power_factor():
    """Plot the thermoelectric power factor with respect to fermi-level"""
    fig, ax = plt.subplots()
    fermi_levels, power_factor, theorical_power_factor = np.genfromtxt("Data/Power factor.csv", unpack=True, delimiter=',', usecols=(0,1,2), skip_header=1)
    ax.plot(fermi_levels * 1e3 / electron_volt, power_factor * 1e3, '-o', markersize=2, c='royalblue', label="Computed")
    # ax.plot(fermi_levels * 1e3 / electron_volt, theorical_power_factor * 1e3, '-o', markersize=2, c='darkorange', label="Bulk/Analytical")
    ax.set_xlabel('Fermi-level (meV)')
    ax.set_ylabel('Power factor (mW/m.$K^{2}$)')
    ax.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    fig.savefig("Power factor.pdf", format="pdf", bbox_inches="tight")
    plt.close(fig)

def plot_mapping_constant():
    """Plot mapping constant with respect to fermi level"""
    material = get_media_class(cf.media)(cf.temp, fermi_level=cf.media_fermi_level)
    fig, ax = plt.subplots()
    fermi_level, mapping_constant = np.genfromtxt("Data/Mapping constant.csv", unpack=True, delimiter=',', usecols=(0,1), skip_header=1)
    material_constant = interpolate_property(fermi_level, mapping_constant, material.fermi_level) * np.ones_like(fermi_level)
    ax.plot(fermi_level * 1e3 / electron_volt, mapping_constant, '-o', markersize=2, c='royalblue', label='Computed')
    ax.plot(fermi_level * 1e3 / electron_volt, material_constant, markersize=2, c='darkorange', label=f"C={material_constant[0]:.4e}")
    ax.set_xlabel('Fermi-level (meV)')
    ax.set_ylabel('Mapping constant (m^2)')
    ax.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    fig.savefig("Mapping constant.pdf", format="pdf", bbox_inches="tight")
    plt.close(fig)

def plot_electron_thermal_conductivity():
    """Plot electron thermal conductivity with respect to fermi-level"""
    fig, ax = plt.subplots()
    fermi_level, thermal_conductivity, theorical_thermal_conductivity = np.genfromtxt("Data/Electron thermal conductivity.csv", unpack=True, delimiter=',', usecols=(0,1,2), skip_header=1)
    ax.plot(fermi_level * 1e3 / electron_volt, thermal_conductivity, '-o', markersize=2, c='royalblue', label="Computed")
    # ax.plot(fermi_level * 1e3 / electron_volt, theorical_thermal_conductivity, '-o', markersize=2, c='darkorange', label="Bulk/Analytical")
    ax.set_xlabel('Fermi-level (meV)')
    ax.set_ylabel('Thermal conductivity (W/m.K)')
    ax.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    fig.savefig("Electron thermal conductivity.pdf", format="pdf", bbox_inches="tight")
    plt.close(fig)

def plot_time_in_segments():
    """Plot time spent in segments"""
    fig, ax = plt.subplots()
    segment, time = np.genfromtxt("Data/Time spent in segments.csv", unpack=True, delimiter=',', usecols=(0, 1), skip_header=1)
    ax.plot(segment, time, '-o', markersize=2, c='royalblue')
    ax.set_xlabel('Y (μm)')
    ax.set_ylabel('Time spent (ns)')
    fig.savefig("Time spent in segments.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)


def plot_thermal_conductivity():
    """Plot thermal conductivity against time segment"""
    fig, ax = plt.subplots()
    time, eff_tc, mat_tc = np.genfromtxt("Data/Thermal conductivity.csv", unpack=True, delimiter=',', usecols=(0, 1, 2), skip_header=1)
    av_eff_tc, av_mat_tc, std_eff_tc, std_mat_tc, av_start, av_end = np.genfromtxt("Data/Average thermal conductivity.csv", unpack=True, delimiter=',', skip_header=1)

    # Plot thermal conductivity vs time:
    ax.plot(time, mat_tc, '-o', markersize=2, linewidth=1, color='deeppink', label='κ$_{mat}$')
    ax.plot(time, eff_tc, '-o', markersize=2, linewidth=1, color='royalblue', label='κ$_{eff}$')

    # Add labels:
    ax.text(av_start/2, max(mat_tc), 'Stabilization period', color='grey', fontsize=8, ha='center')
    ax.text(av_start + (av_end-av_start)/2, max(mat_tc), 'Measurement period', color='grey', fontsize=8, ha='center')

    # Plot the range of averaging:
    ax.axvline(x=av_start, color='grey', linestyle='--')
    ax.axvline(x=av_end, color='grey', linestyle='--')

    # Plot the standard deviation range:
    rectangle_mat = patches.Rectangle((av_start, av_mat_tc-std_mat_tc), av_end-av_start, 2*std_mat_tc, linewidth=0, facecolor='deeppink', alpha=0.2)
    rectangle_eff = patches.Rectangle((av_start, av_eff_tc-std_eff_tc), av_end-av_start, 2*std_eff_tc, linewidth=0, facecolor='royalblue', alpha=0.2)
    ax.add_patch(rectangle_eff)
    ax.add_patch(rectangle_mat)

    # Plot the averaged values in the averaging range:
    ax.plot([av_start, av_end], [av_mat_tc, av_mat_tc], '-', markersize=2, linewidth=2, color='deeppink', label=r'$\overline{\kappa}_{mat}$')
    ax.plot([av_start, av_end], [av_eff_tc, av_eff_tc], '-', markersize=2, linewidth=2, color='royalblue', label=r'$\overline{\kappa}_{eff}$')

    ax.set_ylabel('Thermal conductivity (W/m·K)')
    ax.set_xlabel('Time (ns)')

    ax.set_xlim(left=0.0)
    ax.legend()
    kappa_mat_str = r"$\overline{\kappa}_{mat}$"
    kappa_eff_str = r"$\overline{\kappa}_{eff}$"
    ax.set_title(f'{kappa_mat_str} = {av_mat_tc} ± {std_mat_tc} W/m·K, {kappa_eff_str} = {av_eff_tc} ± {std_eff_tc} W/m·K', color='grey')
    fig.savefig("Thermal conductivity.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)


def plot_temperature_profile():
    """Plot profile of temperature for each time segment"""
    fig, ax = plt.subplots()
    data = np.genfromtxt("Data/Temperature profiles y.csv", unpack=True, delimiter=',', skip_header=1, encoding='utf-8')
    for timeframe_num in range(len(data) - 1):
        ax.plot(data[0][1:], data[timeframe_num + 1][1:], linewidth=1, label=f'Time frame {timeframe_num+1}')
    ax.set_xlabel('Y (μm)')
    ax.set_ylabel('Temperature (K)')
    ax.legend()
    fig.savefig("Temperature profile.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)

def plot_heat_flux_profile():
    """Plot profile of heat flux for each time segment"""

    # Effective heat flux:
    fig, ax = plt.subplots()
    data = np.genfromtxt("Data/Heat flux profiles y.csv", unpack=True, delimiter=',', skip_header=1, encoding='utf-8')
    for timeframe_num in range((len(data) - 1) // 2):
        ax.plot(data[0], data[timeframe_num + 1], linewidth=1, label=f'Time frame {timeframe_num+1}')
    ax.set_xlabel('Y (μm)')
    ax.set_ylabel('Heat flux (W/m²)')
    ax.legend()
    fig.savefig("Heat flux profile effective.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)

    # Material heat flux:
    fig, ax = plt.subplots()
    for timeframe_num in range((len(data) - 1) // 2):
        ax.plot(data[0], data[2*timeframe_num + 1], linewidth=1, label=f'Time frame {timeframe_num+1}')
    ax.set_xlabel('Y (μm)')
    ax.set_ylabel('Heat flux (W/m²)')
    ax.legend()
    fig.savefig("Heat flux profile material.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)

def plot_thermal_map():
    """Plot thermal map as color map"""
    fig = plt.figure()
    thermal_map = np.genfromtxt("Data/Thermal map.csv", unpack=False, delimiter=',', skip_header=0, encoding='utf-8')
    thermal_map = np.flipud(thermal_map)
    minimum_of_colorbar = 1e-20  # Cannot be zero!
    boundaries = [(-cf.width / 2) * 1e6, (cf.width / 2) * 1e6, 0, cf.length * 1e6]
    plt.imshow(thermal_map, cmap='hot', interpolation='none', extent=boundaries,
               norm=LogNorm(vmin=minimum_of_colorbar, vmax=np.amax(thermal_map)))
    plt.xlabel('x (μm)')
    plt.ylabel('y (μm)')
    cbar = plt.colorbar(shrink=0.6)
    cbar.set_label('Energy density', rotation=90)
    fig.savefig("Thermal map.pdf", bbox_inches="tight")
    plt.close(fig)

def plot_pixel_volumes():
    """Plot the pixel volumes as 2D map"""
    fig = plt.figure()
    pixel_volumes = np.genfromtxt("Data/Pixel volumes.csv", unpack=False, delimiter=',', skip_header=0, encoding='utf-8')
    pixel_volumes = np.flipud(pixel_volumes)
    boundaries = [(-cf.width / 2) * 1e6, (cf.width / 2) * 1e6, 0, cf.length * 1e6]
    plt.imshow(pixel_volumes, cmap='hot', interpolation='none', extent=boundaries)
    plt.xlabel('x (μm)')
    plt.ylabel('y (μm)')
    fig.savefig("Pixel volumes.pdf", bbox_inches="tight")
    plt.close(fig)

def plot_heat_flux_map(file, label, units="a.u."):
    """Plot heat flux map as color map"""
    fig = plt.figure()
    heat_flux_map = np.genfromtxt(file, unpack=False, delimiter=',', skip_header=0, encoding='utf-8')
    heat_flux_map = np.flipud(heat_flux_map)
    minimum_of_colorbar = 1e5
    boundaries = [(-cf.width / 2) * 1e6, (cf.width / 2) * 1e6, 0, cf.length * 1e6]
    plt.imshow(heat_flux_map, cmap='jet', interpolation='none', extent=boundaries,
               norm=LogNorm(vmin=minimum_of_colorbar, vmax=np.amax(heat_flux_map)))
    plt.xlabel('x (μm)', fontsize=10)
    plt.ylabel('y (μm)', fontsize=10)
    cbar = plt.colorbar(shrink=0.6)
    cbar.set_label(f"{label} ({units})", rotation=90, fontsize=10)
    fig.savefig(f"{label}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_scattering_map():
    """Plot the map of scattering events"""
    if not cf.output_scattering_map:
        return
    fig, ax = plt.subplots()
    filename = "Data/Scattering map.csv"
    spec_x, spec_y, diff_x, diff_y, int_x, int_y = np.genfromtxt(filename, unpack=True, delimiter=',', skip_header=1, encoding='utf-8')
    ax.plot(np.trim_zeros(diff_x, 'b')*1e6, np.trim_zeros(diff_y, 'b')*1e6, 'o', color='b', markersize=0.1, alpha=0.3)
    ax.plot(np.trim_zeros(spec_x, 'b')*1e6, np.trim_zeros(spec_y, 'b')*1e6, 'o', color='g', markersize=0.1, alpha=0.3)
    ax.plot(np.trim_zeros(int_x, 'b')*1e6, np.trim_zeros(int_y, 'b')*1e6, 'o', color='r', markersize=0.1, alpha=0.3)
    ax.set_xlabel('x (μm)')
    ax.set_ylabel('y (μm)')
    ax.legend(["Diffuse", "Specular", "Internal"])
    ax.set_aspect('equal', 'datalim')
    fig.savefig("Scattering map.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_trajectories():
    """Plot the particle trajectories"""

    data = np.genfromtxt("Data/Particle paths.csv", unpack=False, delimiter=',', skip_header=1, encoding='utf-8')

    # Create XY plot:
    fig, ax = plt.subplots()

    # Draw structure:
    patches = draw_structure_top_view(cf, color_holes='black', color_back=cf.output_structure_color)
    for patch in patches:
        ax.add_patch(patch)


    # Draw paths:
    for index in range(cf.output_trajectories_of_first):
        x_coordinates = np.trim_zeros(data[:, 3 * index], trim='b')
        y_coordinates = np.trim_zeros(data[:, 3 * index + 1], trim='b')
        max_steps = min([x_coordinates.shape[0], y_coordinates.shape[0]])
        ax.plot(x_coordinates[:max_steps], y_coordinates[:max_steps], linewidth=0.2)

    # Set labels:
    ax.set_xlabel('X (μm)')
    ax.set_ylabel('Y (μm)')
    # ax.set_aspect('equal', 'datalim') #remove this to adapt x and y limit 
    ax.set_xlim(2*-cf.width*1e6, 2*cf.width*1e6)
    ax.set_ylim(0, cf.length*1e6)
    fig.savefig("Particle paths XY.pdf", dpi=600, format='pdf', bbox_inches="tight")
    plt.close(fig)

    # Create YZ plot:
    fig, ax = plt.subplots()

    # Draw structure:
    patches = draw_structure_side_view(cf, color_holes='white', color_back=cf.output_structure_color)
    for patch in patches:
        ax.add_patch(patch)

    for index in range(cf.output_trajectories_of_first):
        y_coordinates = np.trim_zeros(data[:, 3 * index + 1], trim='b')
        num_of_points = len(y_coordinates)
        z_coordinates = data[:num_of_points, 3 * index + 2]
        max_steps = min([y_coordinates.shape[0], z_coordinates.shape[0]])
        ax.plot(y_coordinates[:max_steps], z_coordinates[:max_steps], linewidth=0.2)
    ax.set_xlabel('Y (μm)')
    ax.set_ylabel('Z (μm)')
    ax.set_aspect('equal', 'datalim')
    fig.savefig("Particle paths YZ.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)


def plot_scattering_statistics():
    """Calculate and plot rates of different scattering events in each length segment"""
    # Load the data from files:
    filename = "Data/Scattering events statistics.csv"
    scattering_data = np.genfromtxt(filename, unpack=True, delimiter=',', skip_header=1, encoding='utf-8')
    filename = "Data/Time spent in segments.csv"
    segments, time_spent = np.genfromtxt(filename, unpack=True, delimiter=',', skip_header=1, encoding='utf-8')

    # Create the plot:
    fig, ax = plt.subplots()
    all_scattering_rates = []
    for scattering_type in range(scattering_data.shape[0]):
        if scattering_type == 6: # Skip hot side scattering
            continue
        # Calculate and plot scattering events per second:
        scattering_rate = [events / time for events, time in zip(scattering_data[scattering_type, :], time_spent)]
        all_scattering_rates.append(scattering_rate)
        ax.plot(segments, scattering_rate,  '-o', markersize=2)

    ax.set_xlabel('Y (μm)')
    ax.set_ylabel('Scattering rate (1/ns)')
    legend = ["Sidewalls diffuse", "Sidewalls specular", "Top & bottom diffuse", "Top & bottom specular",
              "Holes diffuse", "Holes specular", "Internal", "Pillars diffuse", "Pillars specular"]
    # ax.legend(legend, loc='upper right')
    ax.legend(legend, loc='lower center', ncol=2, fancybox=True)
    plt.yscale('log')
    ax.set_ylim(bottom=1.0)
    fig.savefig("Scattering rates.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)

    # Save the file:
    filename = "Data/Scattering rates.csv"
    header1 = "Y [um], Sidewalls diffuse [1/ns], Sidewalls specular [1/ns], Top & bottom diffuse [1/ns], "
    header2 = "Top & bottom specular [1/ns], Holes diffuse [1/ns], Holes specular [1/ns], Hot side [1/ns], "
    header3 = "Internal [1/ns], Pillars diffuse [1/ns], Pillars specular [1/ns]"
    header = header1 + header2 + header3
    data = np.vstack((segments, all_scattering_rates)).T
    np.savetxt(filename, data, fmt='%1.2e', delimiter=",", header=header)

def plot_scattering_rate_vs_energy():
    fig, ax = plt.subplots()
    energies, scattering_rate = np.genfromtxt("Data/Scattering rate vs energy.csv", unpack=True, delimiter=',', usecols=(0,1), skip_header=1)

    ax.plot(energies * 1e3 / electron_volt, scattering_rate, '-o', markersize=2, c='royalblue')
    ax.axhline(y=cf.timestep, color='gray', linestyle='--', linewidth=1, label='Timestep')
    ax.set_xlabel('Energy (meV)')
    ax.set_ylabel('Internal scattering rate (s)')
    ax.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    fig.savefig("Scattering rate vs energy.pdf", format="pdf", bbox_inches="tight")
    plt.close(fig)

def plot_structure():
    """Plot the structure with all the elements"""

    # Create XY plot:
    fig, ax = plt.subplots()

    # Draw structures:
    patches = draw_structure_top_view(cf, color_holes='black', color_back='royalblue')
    for patch in patches:
        ax.add_patch(patch)

    # Set labels:
    ax.set_xlabel('X (μm)')
    ax.set_ylabel('Y (μm)')
    ax.set_aspect('equal', 'datalim')
    ax.set_xlim(1.2*-cf.width*1e6, 1.2*cf.width*1e6)
    ax.set_ylim(0, cf.length*1e6)
    fig.savefig("Structure XY.pdf", dpi=600, format='pdf', bbox_inches="tight")
    plt.close(fig)

def plot_material_properties():
    """Plot phonon dispersion and display some other material properties"""

    # Initialize the material:
    if cf.media == "Si":
        material = Si(cf.temp)
    elif cf.media == "Ge ":
        material = Ge(cf.temp)
    elif cf.media == "SiC":
        material = SiC(cf.temp)
    elif cf.media == "Graphite":
        material = Graphite(cf.temp)
    else:
        logging.error(f"Material {cf.media} is not supported")
        sys.exit()

    # Plot phonon dispersion:
    fig, ax = plt.subplots()
    for index in range(material.dispersion.shape[1] - 1):
        ax.plot(material.dispersion[:,0], material.dispersion[:,index + 1] * 1e-12, linewidth=1, label=f'{index}')
    ax.set_xlabel('Wavevector (1/m)')
    ax.set_ylabel('Frequency (THz)')
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0)

    # Add material properties:
    ax.set_title(f'{material.name},  T = {cf.temp} K,  C$_p$ = {material.heat_capacity:.3f} J/kg·K,  ρ = {material.density} kg/m³', color="grey")

    fig.savefig("Material properties.pdf", format='pdf', bbox_inches="tight")
    plt.close(fig)


def plot_data(particle_type: ParticleType, cf, mfp_sampling=False):
    """Create plots of various distributions, maps, profiles, and other quantities"""

    # Defines what will be plotted after the phonon simulations:
    phonon_function_list = [
        plot_structure,
        plot_trajectories,
        plot_angle_distribution,
        plot_free_path_distribution,
        plot_frequency_distribution,
        plot_wavelength_distribution,
        plot_travel_time_distribution,
        plot_mean_free_path_distribution,
        plot_velocity_distribution,
        plot_time_in_segments,
        plot_thermal_conductivity,
        plot_temperature_profile,
        plot_heat_flux_profile,
        plot_thermal_map,
        plot_pixel_volumes,
        plot_scattering_statistics,
        plot_scattering_map,
        plot_material_properties,
    ]

    if cf.holes:
        phonon_function_list.extend([
            plot_scattering_angle_distribution,
            ])

    if cf.interfaces:
        phonon_function_list.extend([
            plot_interfaces_angles_distribution,
            plot_transmission_vs_angle,
            plot_transmission_vs_wavelength,
            plot_transmission_vs_frequency,
            plot_transmission_heatmaps_by_mode,
            ])

    # Defines what will be plotted after the electron simulations:
    electron_function_list = [
        plot_structure,
        plot_trajectories,
        plot_angle_distribution,
        plot_free_path_distribution,
        plot_frequency_distribution,
        plot_wavelength_distribution,
        plot_travel_time_distribution,
        plot_mean_free_path_distribution,
        plot_velocity_distribution,
        plot_energy_distribution,
        plot_travel_time_vs_energy,
        plot_transport_function,
        plot_electron_conductivity,
        plot_seebeck_coefficient,
        plot_power_factor,
        plot_mapping_constant,
        plot_electron_thermal_conductivity,
        plot_scattering_rate_vs_energy,
        plot_time_in_segments,
        # plot_thermal_conductivity,
        # plot_temperature_profile,
        # plot_heat_flux_profile,
        # plot_thermal_map,
        plot_pixel_volumes,
        plot_scattering_statistics,
        plot_scattering_map,
    ]

    if cf.holes:
        electron_function_list.extend([
            plot_scattering_angle_distribution,
            ])

    if particle_type is ParticleType.PHONON:
        function_list = phonon_function_list
    else:
        function_list = electron_function_list

    # Run main functions and handle exceptions:
    for func in function_list:
        try:
            func()
        except Exception as e:
            logging.warning(f"Function {func.__name__} failed: {e}")

    # Run additional functions:
    if particle_type is ParticleType.PHONON:
        plot_cumulative_thermal_conductivity(mfp_sampling)
        plot_heat_flux_map(file="Data/Heat flux map xy.csv", label="Heat flux map", units="W/m²")
        plot_heat_flux_map(file="Data/Heat flux map x.csv", label="Heat flux map x", units="W/m²")
        plot_heat_flux_map(file="Data/Heat flux map y.csv", label="Heat flux map y", units="W/m²")




