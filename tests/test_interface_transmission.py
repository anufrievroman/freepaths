"""
VerticalPlane must transmit a phonon with exactly the SMMM probability alpha_total_2T,
evaluated at the TRUE angle of incidence on the plane, i.e. the angle between the phonon
direction and the plane normal x: cos(theta_i) = |d_x| = |cos(phi) sin(theta)|, in 3D.
The fraction of transmitted phonons over many encounters is compared with that probability,
in-plane (phi = 0) and out of the plane, for modes below and above the critical angle.
"""

import math
import random
from types import SimpleNamespace

import numpy as np
import pytest

import freepaths.scatterers as H
from freepaths.interface_smmm import alpha_total_2T
from freepaths.materials import Si, SiGe
from freepaths.options import ParticleType
from freepaths.scattering_types import ScatteringTypes

from surface_cases import make_cf


N = 6000
ROUGHNESS = 0.5e-9


@pytest.fixture(scope="module")
def plane():
    interface = H.VerticalPlane(position_x=0.0, outer_material="Si", inner_material="SiGe")
    interface.assign_materials(SimpleNamespace(temp=300))
    return interface


def expected_transmission(plane, branch, f, d):
    si, sige = plane.outer_material, plane.inner_material
    theta_i = math.acos(abs(d[0]))
    omega = 2 * math.pi * f
    return alpha_total_2T(theta_i, si.group_velocity(branch, f), sige.group_velocity(branch, f),
                          si.density, sige.density, omega, omega, si, sige, branch, ROUGHNESS)


# (branch, f [Hz], in-plane angle from the normal [deg], phi [deg]):
ENCOUNTERS = [
    (0, 1.0e12, 0, 0), (0, 1.0e12, 40, 0), (0, 1.0e12, 60, 0), (1, 1.5e12, 30, 0),
    (0, 1.0e12, 0, 45), (0, 1.0e12, 40, 60), (1, 1.5e12, 30, 50), (0, 3.0e12, 20, 70),
]


@pytest.mark.parametrize("branch, f, in_plane_deg, phi_deg", ENCOUNTERS,
                         ids=[f"branch{b}-{f/1e12:g}THz-{a}deg-phi{p}" for b, f, a, p in ENCOUNTERS])
def test_transmitted_fraction_is_the_smmm_probability(plane, branch, f, in_plane_deg, phi_deg):
    cf = make_cf(interface_roughness=ROUGHNESS)
    theta = math.radians(90 - in_plane_deg)        # theta from the y axis; 90 deg = along +x, the normal
    phi = math.radians(phi_deg)
    d = np.array([math.sin(theta) * math.cos(phi), math.cos(theta) * math.cos(phi), math.sin(phi)])
    expected = expected_transmission(plane, branch, f, d)

    si, sige = plane.outer_material, plane.inner_material
    random.seed(31)
    transmitted = 0
    for _ in range(N):
        pt = SimpleNamespace(type=ParticleType.PHONON, x=-1e-12, y=1e-6, z=0.0, theta=theta, phi=phi,
                             f=f, k=si.wavevector(branch, f), branch_number=branch,
                             speed=si.group_velocity(branch, f), wavelength=2 * math.pi / si.wavevector(branch, f),
                             speed_in=lambda material, b=branch, f=f: material.group_velocity(b, f))
        st = ScatteringTypes()
        plane.scatter(pt, cf, st)
        transmitted += st.interfaces_transmission is not None
    sigma = max(math.sqrt(expected * (1 - expected) / N), 1 / N)
    assert abs(transmitted / N - expected) < 5 * sigma, \
        (f"transmitted {transmitted / N:.4f}, SMMM at the true incidence angle "
         f"{math.degrees(math.acos(abs(d[0]))):.1f} deg gives {expected:.4f}")
