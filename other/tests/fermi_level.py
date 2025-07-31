import numpy as np
from scipy import integrate, optimize, constants, special

# Constants
k_B = constants.k           # Boltzmann constant [J/K]
e = constants.e             # Elementary charge [C]
h = constants.h             # Planck constant [J·s]
m0 = constants.m_e          # Electron rest mass [kg]
k_B_eV = k_B / e            # Boltzmann constant [eV/K]

def Nc(T, m_eff_e):
    """Effective density of states in conduction band [cm^-3]."""
    return 2 * ((2 * np.pi * m_eff_e * m0 * k_B * T) / (h ** 2)) ** 1.5 / 1e6

def Nv(T, m_eff_h):
    """Effective density of states in valence band [cm^-3]."""
    return 2 * ((2 * np.pi * m_eff_h * m0 * k_B * T) / (h ** 2)) ** 1.5 / 1e6

def fermi_dirac_half(eta):
    """Fermi–Dirac integral of order 1/2, numerically stable."""
    gamma_3_2 = np.sqrt(np.pi) / 2
    # on utilise expit(eta - x) pour éviter l'overflow
    integrand = lambda x: np.sqrt(x) * special.expit(eta - x)
    result, _ = integrate.quad(integrand, 0, np.inf)
    return result / gamma_3_2

def fermi_level(n=None, p=None, T=300, Eg=1.12,
                m_eff_e=1.18, m_eff_h=0.81,
                method="fd"):
    """
    Calcule la position du niveau de Fermi :
      - si n fourni    --> renvoie E_C - E_F [eV]
      - si p fourni    --> renvoie E_F - E_V [eV]
    method : "fd" (dégénéré, Fermi-Dirac) ou "boltzmann" (non-dégénéré)
    """
    kT = k_B_eV * T

    # Densités d'états
    Nc_val = Nc(T, m_eff_e)
    Nv_val = Nv(T, m_eff_h)

    if method.lower() == "boltzmann":
        # Approximation non-dégénérée : n = Nc * exp((E_F - E_C)/kT)
        #                              p = Nv * exp((E_V - E_F)/kT)
        if n is not None and p is None:
            EC_minus_EF = kT * np.log(Nc_val / n)
            return EC_minus_EF
        elif p is not None and n is None:
            EF_minus_EV = kT * np.log(Nv_val / p)
            return EF_minus_EV
        else:
            raise ValueError("Spécifiez exactement n ou p, pas les deux.")

    elif method.lower() == "fd":
        # Dégénéré : résolution numérique de Nc * F_{1/2}(eta) = n, etc.
        if n is not None and p is None:
            func = lambda eta: Nc_val * fermi_dirac_half(eta) - n
            eta = optimize.brentq(func, -50, 50)
            EC_minus_EF = -eta * kT
            return EC_minus_EF
        elif p is not None and n is None:
            func = lambda eta: Nv_val * fermi_dirac_half(eta) - p
            eta = optimize.brentq(func, -50, 50)
            EF_minus_EV = -eta * kT
            return EF_minus_EV
        else:
            raise ValueError("Spécifiez exactement n ou p, pas les deux.")

    else:
        raise ValueError("method doit être 'fd' ou 'boltzmann'.")

if __name__ == "__main__":
    method = "fd"
    n = 3e20
    p = 3e20
    # Exemples
    # Non-dégénéré (faible dopage)
    print(f"{method}, n = {n} cm^-3 : {fermi_level(n=n, method=method):.03f}")
    # Dégénéré (fort dopage)
    print(f"{method}, p = {p} cm^-3 : {fermi_level(p=p, method=method):.03f}")
