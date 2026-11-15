"""This module implements calulations of interlayer transmission based on smmm formalizm"""

import math
import numpy as np


def get_k_from_omega(omega, material, branch_number):
    """
    Compute k(ω) from the dispersion relation of the material
    for a given polarization (branch_number).

    - omega : phonon angular frequency (rad/s)
    - material : object containing the dispersion table [k, f_LA, f_TA1, f_TA2]
    - branch_number : 0 (LA), 1 (TA1), 2 (TA2)

    Returns:
    - k : wavevector corresponding to ω
    """
    freq_branch = material.dispersion[:, branch_number + 1]  # ex: colonne 1 pour LA
    k_values = material.dispersion[:, 0]
    omega_values = 2 * np.pi * freq_branch
    return np.interp(omega, omega_values, k_values)


def lambda_from_k(k):
    """Compute the wavelength from the wavevector"""
    return 2 * math.pi / k if k != 0 else float('inf')


def alpha_spec(theta_i, vg_i, vg_j, rho_i, rho_j):
    """
    Specular transmission coefficient (AMM) between two materials.

    Ran & Cao, PRB 110, 024302 (2024), Eq. (20):

        alpha = 4*A*B / (A + B)**2,   A = Z_j/Z_i,  B = cos(theta_j)/cos(theta_i)

    which is the standard oblique-incidence acoustic result 4*W_i*W_j/(W_i + W_j)**2
    with the angle-corrected impedances W = Z/cos(theta), i.e. the COSINES IN THE
    DENOMINATOR ARE CROSSED: (Z_i*cos(theta_j) + Z_j*cos(theta_i))**2. Writing them
    uncrossed agrees at normal incidence but over-estimates transmission at oblique
    incidence (Si/SiGe LA at 1 THz: 0.990 vs 0.928 at 60 deg).

    The transmitted angle follows Snell's law, Eq. (18): vg_i/sin(theta_i) = vg_j/sin(theta_j),
    i.e. sin(theta_j) = (vg_j/vg_i)*sin(theta_i). A phonon entering a slower medium bends
    towards the normal; total reflection happens only on the way into a faster medium, past
    arcsin(vg_i/vg_j).

    Note this is Eq. (20), not Eq. (25): the MTM detailed-balance prefactor
    min{1, P_j*C_j*vg_j**3 / (P_i*C_i*vg_i**3)} of Eq. (25) is deliberately not applied.

    Parameters:
    theta_i -- incidence angle
    vg_i, vg_j -- group velocities of media i and j
    rho_i, rho_j -- mass densities

    Returns:
    Specular transmission coefficient, between 0 and 1
    """
    if vg_i <= 0 or vg_j <= 0:
        return 0.0
    sin_theta_j = (vg_j / vg_i) * math.sin(theta_i)
    if abs(sin_theta_j) >= 1:
        return 0.0  # Total reflection

    theta_j = math.asin(sin_theta_j)
    Z_i = rho_i * vg_i
    Z_j = rho_j * vg_j

    num = 4 * Z_i * Z_j * abs(math.cos(theta_i)) * abs(math.cos(theta_j))
    den = (Z_i * abs(math.cos(theta_j)) + Z_j * abs(math.cos(theta_i))) ** 2
    return num / den


def alpha_diff(k_i, k_j, P_i, P_j):
    """
    Diffuse transmission coefficient based on wavevectors and specularity.

    Parameters:
    k_i, k_j -- incident and transmitted wavevectors
    P_i, P_j -- specularity on side i and j

    Returns:
    Diffuse transmission coefficient, between 0 and 1
    """
    num = (1 - P_j) * k_j**2
    den = (1 - P_i) * k_i**2 + (1 - P_j) * k_j**2
    return num / den if den != 0 else 0.0


def alpha_total_2T(theta_i, vg_i, vg_j, rho_i, rho_j, omega_i, omega_j, material_i, material_j, branch_number, roughness):
    """
    Transmission model for 2 transmission : ex: Si → Ge → Si for layers inside single crystal Si
    Total hybrid transmission coefficient:
    T_total = T1 * T2 with:
    T1 = P1 * a_s1 + (1 - P1) * a_d1   (Si → Ge)
    T2 = P2 * a_s2 + (1 - P2) * a_d2   (Ge → Si)
    """
    # A mode above the top of its branch in the layer has no state to go to:
    if omega_j > 2 * math.pi * material_j.dispersion[:, branch_number + 1].max():
        return 0.0
    if vg_i <= 0 or vg_j <= 0:
        return 0.0

    # First transition: from outer crystal to inside the layer
    k_i = get_k_from_omega(omega_i, material_i, branch_number)
    k_j = get_k_from_omega(omega_j, material_j, branch_number)

    wavelength_i = lambda_from_k(k_i)
    wavelength_j = lambda_from_k(k_j)

    P_i = math.exp(-16 * math.pi**2 * roughness**2 / wavelength_i**2)
    P_j = math.exp(-16 * math.pi**2 * roughness**2 / wavelength_j**2)

    # Snell's law, sin(theta_t) = (vg_j/vg_i)*sin(theta_i): total reflection at the entry is
    # possible only if the layer is the faster medium. At the exit the phonon is back at its
    # incidence angle theta_i, so it is never totally reflected there:
    sin_theta_t1 = (vg_j / vg_i) * math.sin(theta_i)
    if abs(sin_theta_t1) >= 1:
        return 0.0
    theta_t1 = math.asin(sin_theta_t1)

    a_s1 = alpha_spec(theta_i, vg_i, vg_j, rho_i, rho_j)
    a_d1 = alpha_diff(k_i, k_j, P_i, P_j)
    T1 = P_i * a_s1 + (1 - P_i) * a_d1

    # Second transition: From inside the layer to outer crystal
    a_s2 = alpha_spec(theta_t1, vg_j, vg_i, rho_j, rho_i)
    a_d2 = alpha_diff(k_j, k_i, P_j, P_i)
    T2 = P_j * a_s2 + (1 - P_j) * a_d2

    return T1 * T2
