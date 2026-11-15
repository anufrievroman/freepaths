"""
End-to-end benchmarks: the real freepaths CLI is run on small input files in a temporary
folder, and its outputs are compared with exact limits. Slow (minutes); select or deselect
them with `-m slow` / `-m "not slow"`.

1. Diffusive limit: Fourier tracing in a "gray" material (one group velocity v, one relaxation
   time tau, see gray_material.py), in a box 20 mean free paths long with specular walls, must
   give the bulk kappa = C v^2 tau / 3. The only analytic check of the whole Fourier pipeline:
   emission, rethermalization, energy and flux maps, the heat capacity, and the gradient fit.
2. Ballistic limit: gray material, no internal scattering, specular walls. The interior then sits
   at half the reservoir temperature difference, so J / (2 T_interior) must equal the Landauer
   conductance (1/4) sum D C v of Material.ballistic_conductance(). Checks the emission spectrum,
   the Lambert injection, and the flux and temperature bookkeeping together.
3. Clausing transmission: a duct with fully diffuse walls, no internal scattering, absorbing
   ends. The fraction of particles from the hot end that reach the cold end must match an
   independent ray tracer, which itself reproduces Clausing's published values for a
   cylindrical tube (FreePATHS cannot build a cylinder).
"""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import gray_material
import two_tau_material


pytestmark = pytest.mark.slow

REPO = Path(__file__).resolve().parent.parent
TESTS = Path(__file__).resolve().parent


def run_freepaths(folder, name, body):
    """Run `freepaths <input>` in the folder with the working tree on PYTHONPATH; returns the Data folder"""
    input_file = folder / f"{name}.py"
    input_file.write_text(f"OUTPUT_FOLDER_NAME = {name!r}\nLOW_MEMORY_USAGE = True\n"
                          f"OUTPUT_TRAJECTORIES_OF_FIRST = 0\nNUMBER_OF_PROCESSES = 8\n" + body)
    env = dict(os.environ, PYTHONPATH=f"{REPO}{os.pathsep}{TESTS}")
    result = subprocess.run([sys.executable, "-m", "freepaths", str(input_file)], cwd=folder, env=env,
                            capture_output=True, text=True, timeout=3600)
    assert result.returncode == 0, f"freepaths failed:\n{result.stdout[-2000:]}\n{result.stderr[-2000:]}"
    return folder / "Results" / name / "Data"


def interior_profiles(data, frames, stabilization, window=(0.2, 0.8)):
    """Mean temperature and effective heat flux over the interior pixels and the steady timeframes"""
    temperature = np.loadtxt(data / "Temperature profiles y.csv", delimiter=",")
    flux = np.loadtxt(data / "Heat flux profiles y.csv", delimiter=",")
    n = len(temperature)
    rows = slice(int(window[0] * n), int(window[1] * n))
    steady = slice(1 + stabilization, 1 + frames)
    return temperature[rows, steady], flux[rows, steady]


def test_diffusive_limit_of_a_gray_medium(tmp_path):
    box = """
import sys
sys.path.insert(0, {tests!r})
from gray_material import install
install()
MEDIA = "Gray"
T = 300
NUMBER_OF_PARTICLES = 6000
TIMESTEP = 1e-12
NUMBER_OF_TIMESTEPS = 40000
NUMBER_OF_VIRTUAL_TIMESTEPS = 80000
NUMBER_OF_TIMEFRAMES = 8
NUMBER_OF_STABILIZATION_TIMEFRAMES = 4
WIDTH = 300e-9
THICKNESS = 300e-9
LENGTH = 1200e-9
NUMBER_OF_PIXELS_X = 1
NUMBER_OF_PIXELS_Y = 40
GRADIENT_FIT_RANGE = (0.2, 0.8)
INCLUDE_INTERNAL_SCATTERING = True
SIDE_WALL_ROUGHNESS = 0.0
TOP_ROUGHNESS = 0.0
BOTTOM_ROUGHNESS = 0.0
PARTICLE_SOURCES = [Source(x=0, y=0, z=0, size_x=WIDTH, size_y=0, size_z=THICKNESS, angle_distribution="lambert", angle=0)]
""".format(tests=str(TESTS))
    data = run_freepaths(tmp_path, "gray_diffusive", box)
    kappa, _, error, _, _, _ = np.loadtxt(data / "Average thermal conductivity.csv", delimiter=",")
    expected = gray_material.bulk_conductivity(300)
    assert kappa == pytest.approx(expected, rel=0.05), \
        f"Fourier kappa {kappa:.1f} +/- {error:.1f} W/mK, gray bulk C v^2 tau / 3 = {expected:.1f} W/mK"


def test_diffusive_limit_of_a_medium_without_mode_mixing(tmp_path):
    """Two relaxation times, all scattering elastic (two_tau_material.py): the modes never mix, so
    kappa = sum C_i v^2 tau_i / 3 only if the hot side absorbs the phonons that come back and emits
    new ones from the emission distribution. Re-emitting them in the same mode gives the
    harmonic-type average instead, 37% lower here, which no gray medium can show."""
    box = """
import sys
sys.path.insert(0, {tests!r})
from two_tau_material import install
install()
MEDIA = "TwoTau"
T = 300
NUMBER_OF_PARTICLES = 4000
TIMESTEP = 0.5e-12
NUMBER_OF_TIMESTEPS = 300000
NUMBER_OF_VIRTUAL_TIMESTEPS = 600000
NUMBER_OF_TIMEFRAMES = 8
NUMBER_OF_STABILIZATION_TIMEFRAMES = 4
WIDTH = 300e-9
THICKNESS = 300e-9
LENGTH = 1200e-9
NUMBER_OF_PIXELS_X = 1
NUMBER_OF_PIXELS_Y = 40
GRADIENT_FIT_RANGE = (0.2, 0.8)
INCLUDE_INTERNAL_SCATTERING = True
SIDE_WALL_ROUGHNESS = 0.0
TOP_ROUGHNESS = 0.0
BOTTOM_ROUGHNESS = 0.0
PARTICLE_SOURCES = [Source(x=0, y=0, z=0, size_x=WIDTH, size_y=0, size_z=THICKNESS, angle_distribution="lambert", angle=0)]
""".format(tests=str(TESTS))
    data = run_freepaths(tmp_path, "two_tau_diffusive", box)
    kappa, _, error, _, _, _ = np.loadtxt(data / "Average thermal conductivity.csv", delimiter=",")
    expected = two_tau_material.bulk_conductivity(300)
    harmonic = two_tau_material.harmonic_conductivity(300)
    assert kappa == pytest.approx(expected, rel=0.10), \
        (f"Fourier kappa {kappa:.1f} +/- {error:.1f} W/mK, sum C v^2 tau / 3 = {expected:.1f} W/mK "
         f"(a hot side that keeps the mode of the returning phonons gives {harmonic:.1f})")


def test_ballistic_limit_is_the_landauer_conductance(tmp_path):
    """In the gray material (one speed), so that no slow zone-edge modes make the energy density a
    heavy-tailed estimate: each particle deposits energy in proportion to its residence time
    L / (v cos theta), and with real Si dispersions the rare slow ones dominate the mean."""
    frames, stabilization = 8, 4
    box = f"""
import sys
sys.path.insert(0, {str(TESTS)!r})
from gray_material import install
install()
MEDIA = "Gray"
T = 300
NUMBER_OF_PARTICLES = 20000
TIMESTEP = 1e-12
NUMBER_OF_TIMESTEPS = 20000
NUMBER_OF_VIRTUAL_TIMESTEPS = 40000
NUMBER_OF_TIMEFRAMES = {frames}
NUMBER_OF_STABILIZATION_TIMEFRAMES = {stabilization}
WIDTH = 300e-9
THICKNESS = 300e-9
LENGTH = 300e-9
NUMBER_OF_PIXELS_X = 1
NUMBER_OF_PIXELS_Y = 30
INCLUDE_INTERNAL_SCATTERING = False
SIDE_WALL_ROUGHNESS = 0.0
TOP_ROUGHNESS = 0.0
BOTTOM_ROUGHNESS = 0.0
PARTICLE_SOURCES = [Source(x=0, y=0, z=0, size_x=WIDTH, size_y=0, size_z=THICKNESS, angle_distribution="lambert", angle=0)]
"""
    data = run_freepaths(tmp_path, "ballistic", box)
    temperature, flux = interior_profiles(data, frames, stabilization)
    profile = temperature.mean(axis=1)            # time-averaged over the steady timeframes
    assert np.std(profile) < 0.03 * np.mean(profile), "a ballistic interior must have a flat temperature"
    measured = np.mean(flux) / (2 * np.mean(temperature))
    expected = gray_material.Gray(300).ballistic_conductance()
    assert measured == pytest.approx(expected, rel=0.03), \
        f"J / (2 T_interior) = {measured:.4e} W/m^2K, Landauer conductance = {expected:.4e} W/m^2K"


# ---------------------------------------------------------------- boundary-limited conduction

GRAY_BOX = """
import sys
sys.path.insert(0, {tests!r})
from gray_material import install
install()
MEDIA = "Gray"
T = 300
NUMBER_OF_PARTICLES = {particles}
TIMESTEP = 0.5e-12
NUMBER_OF_TIMESTEPS = 200000
NUMBER_OF_VIRTUAL_TIMESTEPS = 400000
NUMBER_OF_TIMEFRAMES = 8
NUMBER_OF_STABILIZATION_TIMEFRAMES = 4
WIDTH = {width}
THICKNESS = {thickness}
LENGTH = {length}
NUMBER_OF_PIXELS_X = 1
NUMBER_OF_PIXELS_Y = 40
GRADIENT_FIT_RANGE = (0.2, 0.8)
INCLUDE_INTERNAL_SCATTERING = {internal}
SIDE_WALL_ROUGHNESS = {side_roughness}
TOP_ROUGHNESS = 1e-6
BOTTOM_ROUGHNESS = 1e-6
PARTICLE_SOURCES = [Source(x=0, y=0, z=0, size_x=WIDTH, size_y=0, size_z=THICKNESS, angle_distribution="lambert", angle=0)]
"""

ROUGH = 1e-6    # [m] roughness that makes every wall hit diffuse (Soffer p ~ 0 for all gray modes)


def test_casimir_limit_of_a_diffuse_square_wire(tmp_path):
    """No internal scattering, fully diffuse walls: kappa = C v Lambda / 3 with the Casimir length
    Lambda = 1.115 a of a square wire of side a (McCurdy, Maris & Elbaum, Phys. Rev. B 2, 4077
    (1970)). It is NOT the mean chord 4V/S = a: long flights along the axis carry more heat."""
    side = 50e-9
    box = GRAY_BOX.format(tests=str(TESTS), particles=6000, width=side, thickness=side, length=2000e-9,
                          internal=False, side_roughness=ROUGH)
    data = run_freepaths(tmp_path, "casimir_wire", box)
    kappa, _, error, _, _, _ = np.loadtxt(data / "Average thermal conductivity.csv", delimiter=",")
    expected = gray_material.Gray(300).dispersion_heat_capacity * gray_material.VELOCITY * 1.115 * side / 3
    assert kappa == pytest.approx(expected, rel=0.07), \
        f"Fourier kappa {kappa:.2f} +/- {error:.2f} W/mK, Casimir C v (1.115 a) / 3 = {expected:.2f} W/mK"


def fuchs_sondheimer(delta):
    """kappa_film / kappa_bulk of a film of thickness delta mean free paths with diffuse surfaces:
    1 - (3 / 2 delta) int_1^inf (1/x^3 - 1/x^5) (1 - exp(-delta x)) dx  (Sondheimer, Adv. Phys. 1, 1 (1952))"""
    from scipy.integrate import quad
    integral, _ = quad(lambda x: (x**-3 - x**-5) * (1 - np.exp(-delta * x)), 1, np.inf)
    return 1 - 3 / (2 * delta) * integral


def test_fuchs_sondheimer_reduction_in_a_thin_film(tmp_path):
    """Gray medium in a film one mean free path thick, diffuse top and bottom, specular side walls:
    the in-plane kappa must be the bulk value times the Fuchs-Sondheimer factor."""
    assert fuchs_sondheimer(1e3) == pytest.approx(1 - 3 / 8e3, rel=1e-4)      # thick-film limit 1 - 3/(8 delta)
    mean_free_path = gray_material.VELOCITY * gray_material.RELAXATION_TIME
    box = GRAY_BOX.format(tests=str(TESTS), particles=6000, width=300e-9, thickness=mean_free_path, length=1200e-9,
                          internal=True, side_roughness=0.0)
    data = run_freepaths(tmp_path, "thin_film", box)
    kappa, _, error, _, _, _ = np.loadtxt(data / "Average thermal conductivity.csv", delimiter=",")
    expected = gray_material.bulk_conductivity(300) * fuchs_sondheimer(1.0)
    assert kappa == pytest.approx(expected, rel=0.06), \
        (f"Fourier kappa {kappa:.1f} +/- {error:.1f} W/mK, Fuchs-Sondheimer {fuchs_sondheimer(1.0):.3f} x bulk "
         f"= {expected:.1f} W/mK")


# ---------------------------------------------------------------- Clausing transmission

def lambert_directions(n, rng, normal_axis=1):
    mu = np.sqrt(rng.random(n))
    psi = 2 * np.pi * rng.random(n)
    s = np.sqrt(1 - mu**2)
    d = np.empty((n, 3))
    others = [a for a in range(3) if a != normal_axis]
    d[:, normal_axis] = mu
    d[:, others[0]] = s * np.cos(psi)
    d[:, others[1]] = s * np.sin(psi)
    return d


def reference_transmission(length, cross_section, size, n=200000, seed=31):
    """Independent ray tracer for a duct along y with fully diffuse (Lambert) walls and absorbing ends.
    cross_section "square" (side = size) or "circle" (radius = size). Returns the transmitted fraction."""
    rng = np.random.default_rng(seed)
    if cross_section == "square":
        p = np.column_stack((rng.uniform(-size/2, size/2, n), np.zeros(n), rng.uniform(-size/2, size/2, n)))
    else:
        r = size * np.sqrt(rng.random(n))
        a = 2 * np.pi * rng.random(n)
        p = np.column_stack((r * np.cos(a), np.zeros(n), r * np.sin(a)))
    d = lambert_directions(n, rng)
    transmitted = 0
    while len(p):
        with np.errstate(divide="ignore", invalid="ignore"):
            t_end = np.where(d[:, 1] > 0, (length - p[:, 1]) / d[:, 1], np.where(d[:, 1] < 0, -p[:, 1] / d[:, 1], np.inf))
            if cross_section == "square":
                t_x = np.where(d[:, 0] > 0, (size/2 - p[:, 0]) / d[:, 0], np.where(d[:, 0] < 0, (-size/2 - p[:, 0]) / d[:, 0], np.inf))
                t_z = np.where(d[:, 2] > 0, (size/2 - p[:, 2]) / d[:, 2], np.where(d[:, 2] < 0, (-size/2 - p[:, 2]) / d[:, 2], np.inf))
                t_wall = np.minimum(t_x, t_z)
            else:
                a2 = d[:, 0]**2 + d[:, 2]**2
                b = 2 * (p[:, 0] * d[:, 0] + p[:, 2] * d[:, 2])
                c = p[:, 0]**2 + p[:, 2]**2 - size**2
                t_wall = np.where(a2 > 0, (-b + np.sqrt(np.maximum(b**2 - 4 * a2 * c, 0))) / (2 * a2), np.inf)
        ends = t_end <= t_wall
        transmitted += np.sum(ends & (d[:, 1] > 0))
        keep = ~ends
        p = p[keep] + t_wall[keep, None] * d[keep]
        # New Lambert direction about the inward wall normal:
        if cross_section == "square":
            normal = np.zeros_like(p)
            on_x = np.abs(np.abs(p[:, 0]) - size/2) < np.abs(np.abs(p[:, 2]) - size/2)
            normal[on_x, 0] = -np.sign(p[on_x, 0])
            normal[~on_x, 2] = -np.sign(p[~on_x, 2])
        else:
            normal = np.column_stack((-p[:, 0], np.zeros(len(p)), -p[:, 2])) / size
        t1 = np.column_stack((np.zeros(len(p)), np.ones(len(p)), np.zeros(len(p))))   # y is always tangent
        t2 = np.cross(normal, t1)
        local = lambert_directions(len(p), rng)       # mu along axis 1
        d = local[:, [1]] * normal + local[:, [0]] * t1 + local[:, [2]] * t2
    return transmitted / n


@pytest.mark.parametrize("l_over_r, clausing", [(1.0, 0.67), (2.0, 0.51), (4.0, 0.36)])
def test_reference_tracer_reproduces_clausing(l_over_r, clausing):
    """Clausing (1932) via Table I of Hall et al., arXiv:1302.2963; values given to two decimals"""
    assert reference_transmission(l_over_r, "circle", 1.0) == pytest.approx(clausing, abs=0.01)


@pytest.mark.parametrize("l_over_a", [1.0, 2.0])
def test_clausing_transmission_of_a_diffuse_duct(tmp_path, l_over_a):
    side = 100e-9
    duct = f"""
MEDIA = "Si"
T = 300
NUMBER_OF_PARTICLES = 20000
TIMESTEP = 1e-12
NUMBER_OF_TIMESTEPS = 300000
NUMBER_OF_VIRTUAL_TIMESTEPS = 300000
WIDTH = {side}
THICKNESS = {side}
LENGTH = {l_over_a * side}
NUMBER_OF_PIXELS_X = 1
NUMBER_OF_PIXELS_Y = 10
INCLUDE_INTERNAL_SCATTERING = False
RETHERMALIZATION_ON_HOT_SIDES = False
SIDE_WALL_ROUGHNESS = 1e-6
TOP_ROUGHNESS = 1e-6
BOTTOM_ROUGHNESS = 1e-6
PARTICLE_SOURCES = [Source(x=0, y=0, z=0, size_x=WIDTH, size_y=0, size_z=THICKNESS, angle_distribution="lambert", angle=0)]
"""
    data = run_freepaths(tmp_path, f"clausing_{l_over_a:g}", duct)
    travel_times = np.loadtxt(data / "All travel times.csv")
    measured = np.mean(travel_times > 0)
    expected = reference_transmission(l_over_a, "square", 1.0)
    sigma = np.sqrt(expected * (1 - expected) / len(travel_times))
    assert abs(measured - expected) < 5 * sigma + 0.005, \
        f"L/a = {l_over_a}: {measured:.4f} of the particles reached the cold side, reference {expected:.4f}"
