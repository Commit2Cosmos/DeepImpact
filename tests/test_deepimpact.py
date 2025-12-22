import pandas as pd
import numpy as np
import os

from pytest import fixture, mark

from deepimpact.solver import Planet

# def test_altitude_decrease(planet):
#     """
#     Test that altitude decreases over time
#     """
#     params = {
#         "radius": 60.0,
#         "velocity": 19000.0,
#         "density": 3000.0,
#         "strength": 1e7,
#         "angle": 45.0,
#         "init_altitude": 100e3,
#         "dt": 0.001,
#         "radians": False,
#     }

#     result = planet.solve_atmospheric_entry(**params)

#     # Basic checks
#     assert isinstance(result, pd.DataFrame), "Result should be a pandas DataFrame"
#     assert len(result) > 0, "Result should contain at least one row"

#     # Check that altitude is always non-negative (with tolerance for numerical precision)
#     altitudes = result["altitude"].values
#     assert np.all(
#         altitudes >= -1e-10
#     ), "Altitude should always be non-negative (within numerical precision)"

#     # Check that altitude generally decreases (allow for small numerical errors)
#     altitude_diffs = np.diff(altitudes)

#     # Count how many steps are non-increasing (decreasing or zero)
#     non_increasing_steps = np.sum(
#         altitude_diffs <= 1e-10
#     )  # Allow small numerical increases
#     total_steps = len(altitude_diffs)

#     # At least 90% of steps should be non-increasing
#     # Allow some tolerance for numerical oscillations
#     assert (
#         non_increasing_steps / total_steps >= 0.8
#     ), "Most altitude steps should be non-increasing"


#     print("Altitude decrease test passed")
#     print(f"Final altitude: {altitudes[-1]:.2e}")
#     print(
#         f"Non-increasing steps: {non_increasing_steps}/{total_steps} = {non_increasing_steps/total_steps:.1%}"
#     )


# Test the mass is decreasing
def test_mass_nonincreasing():
    planet = Planet()
    df = planet.solve_atmospheric_entry(
        radius=10, velocity=20000, density=3000, strength=1e6, angle=45
    )

    m = df["mass"].values
    diffs = np.diff(m)

    #  Mass must never be negative
    assert np.all(m >= -1e-12), f"Mass should never be negative, but min(m) = {m.min()}"

    #  Most steps should show decreasing mass
    small_tol = 1e-6  # numerical tolerance
    non_increasing_steps = np.sum(diffs <= small_tol)
    total_steps = len(diffs)

    assert non_increasing_steps / total_steps >= 0.95, (
        f"Mass should be non-increasing for most steps: "
        f"{non_increasing_steps}/{total_steps} steps"
    )

    #  No significant mass increase
    significant_growth = diffs[diffs > 1e-4]
    assert (
        len(significant_growth) == 0
    ), f"Mass increased significantly in some steps: {significant_growth[:10]}"

    # Mass should behave smoothly near impact
    if len(m) > 5:
        tail_diff = np.diff(m[-5:])
        assert np.all(
            tail_diff <= small_tol
        ), "Mass should decrease or stabilize near impact."

    print("Mass non-increasing test passed.")


# the angle change
def test_angle_curve_shape():
    planet = Planet()
    df = planet.solve_atmospheric_entry(
        radius=20, velocity=18000, density=3000, strength=5e6, angle=45
    )

    theta = df["angle"].values
    diffs = np.diff(theta)
    min_idx = np.argmin(theta)

    # Angle must remain within valid physical range
    assert np.all(theta > -1e-6) and np.all(
        theta < 90 + 1e-6
    ), "Angle must remain between 0 and 90 degrees."

    # Before the minimum, angle should generally decrease
    before = diffs[:min_idx]
    frac_decreasing = np.sum(before <= 1e-6) / len(before)
    assert frac_decreasing >= 0.95, "Angle should mostly decrease before the minimum."

    # After the minimum, angle should generally increase
    after = diffs[min_idx:]
    frac_increasing = np.sum(after >= -1e-6) / len(after)
    assert frac_increasing >= 0.95, "Angle should mostly increase after the minimum."

    # No sudden unrealistic jumps
    max_jump = np.max(np.abs(diffs))
    assert (
        max_jump < 2.0
    ), f"Angle changed too abruptly: jump={max_jump} degrees per step."

    # Angle at impact must be close to 90
    if df["altitude"].iloc[-1] <= 1e-6:
        assert (
            theta[-1] > 70
        ), f"Impact angle should be steep, got {theta[-1]:.2f} degrees."


def test_radius_behavior():
    planet = Planet()
    df = planet.solve_atmospheric_entry(
        radius=20, velocity=18000, density=3000, strength=5e6, angle=45
    )

    r = df["radius"].values
    dr = np.diff(r)

    assert (r >= 0).all(), "Radius must be non-negative"
    assert not np.any(dr < -1e-9), "Radius should not decrease"


# spreading_rate
def test_spreading_rate():
    planet = Planet()
    df = planet.solve_atmospheric_entry(
        radius=20, velocity=18000, density=3000, strength=5e6, angle=45
    )

    u = df["spreading_rate"].values

    assert (u >= -1e-9).all(), "Spreading rate should never be negative"


@fixture(scope="module")
def deepimpact():
    import deepimpact

    return deepimpact


@fixture(scope="module")
def planet(deepimpact):
    return deepimpact.Planet()


@fixture(scope="module")
def loc(deepimpact):
    return deepimpact.GeospatialLocator()


@fixture(scope="module")
def result(planet):
    input = {
        "radius": 1.0,
        "velocity": 2.0e4,
        "density": 3000.0,
        "strength": 1e5,
        "angle": 30.0,
        "init_altitude": 100e3,
    }

    result = planet.solve_atmospheric_entry(**input)

    return result


@fixture(scope="module")
def outcome(planet, result):
    energy = planet.calculate_energy(result=result)
    outcome = planet.analyse_outcome(result=energy)
    return outcome


# ============================================================================
# ORIGINAL TESTS
# ============================================================================


def test_import(deepimpact):
    assert deepimpact


def test_planet_signature(deepimpact):
    inputs = {
        "atmos_func": "exponential",
        "atmos_filename": None,
        "Cd": 1.0,
        "Ch": 0.1,
        "Q": 1e7,
        "Cl": 1e-3,
        "alpha": 0.3,
        "Rp": 6371e3,
        "g": 9.81,
        "H": 8000.0,
        "rho0": 1.2,
    }

    # initialise using keyword arguments
    planet_kw = deepimpact.Planet(**inputs)
    # initialise using positional arguments
    planet_pos = deepimpact.Planet(*inputs.values())

    # Check all attributes are set correctly
    for key in inputs.keys():
        if key[:5] == "atmos":
            # you may want to modify this to check that
            # atmos_filename is set correctly
            continue
        else:
            assert planet_kw.__dict__[key] == inputs[key]
            assert planet_pos.__dict__[key] == inputs[key]


def test_attributes(planet):
    for key in ("Cd", "Ch", "Q", "Cl", "alpha", "Rp", "g", "H", "rho0"):
        assert hasattr(planet, key)


def test_atmos_filename(planet):

    assert os.path.isfile(planet.atmos_filename)


def test_solve_atmospheric_entry(result):

    assert type(result) is pd.DataFrame

    for key in ("velocity", "mass", "angle", "altitude", "distance", "radius", "time"):
        assert key in result.columns


def test_calculate_energy(planet, result):

    energy = planet.calculate_energy(result=result)

    assert type(energy) is pd.DataFrame

    for key in (
        "velocity",
        "mass",
        "angle",
        "altitude",
        "distance",
        "radius",
        "time",
        "dedz",
    ):
        assert key in energy.columns


def test_analyse_outcome(outcome):

    assert type(outcome) is dict

    for key in (
        "outcome",
        "burst_peak_dedz",
        "burst_altitude",
        "burst_distance",
        "burst_energy",
    ):
        assert key in outcome.keys()


def test_analytic(planet):
    """
    Tests the numerical solution from solver.py implementation against given data from scenario.npz.
    """
    inputs = {
        "radius": 35.0,
        "angle": 45.0,
        "strength": 1e7,
        "density": 3000.0,
        "velocity": 19e3,
        "init_altitude": 100e3,
    }
    data = np.load("./tests/scenario.npz")
    result2 = pd.DataFrame({x: data[x] for x in data.files})
    result = planet.solve_atmospheric_entry(
        radius=inputs["radius"],
        angle=inputs["angle"],
        strength=inputs["strength"],
        density=inputs["density"],
        velocity=inputs["velocity"],
        init_altitude=inputs["init_altitude"],
        radians=False,
    )
    result = result[result["altitude"] >= result2["altitude"].iloc[-1]]

    n = 40
    indices = np.linspace(0, len(result) - 1, n, dtype=int)

    result = result.iloc[indices]

    result = planet.calculate_energy(result)
    result2 = planet.calculate_energy(result2)

    common_columns = result2.columns
    result = result[common_columns]
    np.testing.assert_allclose(result, result2, rtol=2e-2, atol=1e3)


def test_mass_angle_constraints(planet):
    """This test checks that mass continuously decreases and remains >0 throughout the sim
    We also add a new check that angle continuously increases throughout the sim."""
    inputs = {
        "radius": 35.0,
        "angle": 45.0,
        "strength": 1e7,
        "density": 3000.0,
        "velocity": 19e3,
    }
    result = planet.solve_atmospheric_entry(
        init_altitude=100e3, dt=0.1, radians=False, **inputs
    )

    mass_values = result["mass"].values
    angle_values = result["angle"].values  # Angles are in degrees from the solver

    # Mass not negative
    assert np.all(mass_values >= -1e-9)

    # 2. Mass must always decrease therefore mass only decrease or stay constant mass differences should be <= 0
    mass_diff = np.diff(mass_values)

    assert np.all(mass_diff <= 1e-9)

    # Angle must always increase therefore angle differences should be >= 0
    angle_diff = np.diff(angle_values)

    # Find the minimum angle index
    min_idx = np.argmin(angle_values)

    # Before the minimum: angle should decrease (or stay almost constant)
    assert np.all(
        angle_diff[:min_idx] <= 1e-9
    ), "Angle should decrease before the minimum."

    # After the minimum: angle should increase (or stay almost constant)
    assert np.all(
        angle_diff[min_idx:] >= -1e-9
    ), "Angle should increase after the minimum."

    # Angle must be non-negative

    assert np.all(angle_values >= -1e-9)

    assert np.all(angle_values <= 180.1)


def test_atmos_func_tabular(deepimpact):
    """This test checks that the values provided in the table match
    the ones calculated by the interpolation function."""
    planet_tabular = deepimpact.Planet(atmos_func="tabular")

    alt_rho = pd.read_csv(
        os.sep.join(
            (os.path.dirname(__file__), "..", "resources", "AltitudeDensityTable.csv")
        ),
        sep=r"\s+",
        skiprows=1,
        names=["altitude", "rho"],
    )

    calculated_rho = [planet_tabular.rhoa(x) for x in alt_rho["altitude"]]
    assert np.allclose(alt_rho["rho"].to_numpy(), np.array(calculated_rho))

# ============================================================================
# PARAMETRIZED TESTS FOR INDIVIDUAL PARAMETERS
# ============================================================================


@mark.parametrize("radius", [0.1, 1.0, 5.0, 10.0, 35.0, 50.0, 100.0])
def test_radius_parameter(planet, radius):
    """Test different meteor radii"""
    input = {
        "radius": radius,
        "velocity": 2.0e4,
        "density": 3000.0,
        "strength": 1e5,
        "angle": 30.0,
        "init_altitude": 100e3,
    }
    result = planet.solve_atmospheric_entry(**input)

    # Assertions
    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0
    assert np.isclose(result["radius"].iloc[0], radius, rtol=1e-2)


@mark.parametrize("velocity", [1e2, 1e4, 1.5e4, 2e4, 3e4, 5e4, 7e4])
def test_velocity_parameter(planet, velocity):
    """Test different impact velocities"""
    input = {
        "radius": 10.0,
        "velocity": velocity,
        "density": 3000.0,
        "strength": 1e5,
        "angle": 30.0,
        "init_altitude": 100e3,
    }
    result = planet.solve_atmospheric_entry(**input)

    # Generate plot (uncomment import at top and uncomment line below)
    # plot_parameter_test(result, "velocity", velocity)

    # Assertions
    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0
    assert np.isclose(result["velocity"].iloc[0], velocity, rtol=1e-2)


@mark.parametrize("density", [1000.0, 2000.0, 3000.0, 5000.0, 8000.0])
def test_density_parameter(planet, density):
    """Test different meteor densities (ice to iron)"""
    input = {
        "radius": 10.0,
        "velocity": 2.0e4,
        "density": density,
        "strength": 1e5,
        "angle": 30.0,
        "init_altitude": 100e3,
    }
    result = planet.solve_atmospheric_entry(**input)

    # Assertions
    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0
    initial_mass = (4 / 3) * np.pi * 10.0**3 * density
    assert np.isclose(result["mass"].iloc[0], initial_mass, rtol=1e-2)


@mark.parametrize(
    "mass_target",
    [1e1, 1e3, 1e7, 1e9, 1e11],
)
def test_mass_parameter(planet, mass_target):
    """Test different meteor masses"""
    density = 3000.0
    # Calculate radius to achieve exact target mass
    radius = (3 * mass_target / (4 * np.pi * density)) ** (1 / 3)

    input = {
        "radius": radius,
        "velocity": 2.0e4,
        "density": density,
        "strength": 1e5,
        "angle": 30.0,
        "init_altitude": 100e3,
    }
    result = planet.solve_atmospheric_entry(**input)

    # Assertions
    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0
    assert np.isclose(result["mass"].iloc[0], mass_target, rtol=1e-3)
    # Verify mass conservation or ablation
    assert result["mass"].iloc[-1] <= result["mass"].iloc[0]


@mark.parametrize("strength", [1e4, 1e5, 1e6, 1e7, 1e8])
def test_strength_parameter(planet, strength):
    """Test different material strengths"""
    input = {
        "radius": 10.0,
        "velocity": 2.0e4,
        "density": 3000.0,
        "strength": strength,
        "angle": 30.0,
        "init_altitude": 100e3,
    }
    result = planet.solve_atmospheric_entry(**input)

    # Assertions
    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0


@mark.parametrize("angle", [10.0, 20.0, 30.0, 45.0, 60.0, 75.0, 85.0])
def test_angle_parameter(planet, angle):
    """Test different entry angles"""
    input = {
        "radius": 10.0,
        "velocity": 2.0e4,
        "density": 3000.0,
        "strength": 1e5,
        "angle": angle,
        "init_altitude": 100e3,
    }
    result = planet.solve_atmospheric_entry(**input)

    # Assertions
    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0
    assert np.isclose(result["angle"].iloc[0], angle, rtol=1e-2)


@mark.parametrize("init_altitude", [50e3, 75e3, 100e3, 150e3, 200e3])
def test_init_altitude_parameter(planet, init_altitude):
    """Test different initial altitudes"""
    input = {
        "radius": 10.0,
        "velocity": 2.0e4,
        "density": 3000.0,
        "strength": 1e5,
        "angle": 30.0,
        "init_altitude": init_altitude,
    }
    result = planet.solve_atmospheric_entry(**input)

    # Assertions
    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0
    assert result["altitude"].iloc[0] <= init_altitude


# ============================================================================
# PARAMETRIZED TESTS
# ============================================================================


@mark.parametrize(
    "radius,velocity",
    [
        (0.5, 1e4),  # Small, slow
        (1.0, 2e4),  # Small, medium
        (10.0, 2e4),  # Medium, medium
        (35.0, 3e4),  # Large, fast
        (50.0, 5e4),  # Very large, very fast
    ],
)
def test_radius_velocity_combinations(planet, radius, velocity):
    """Test realistic radius and velocity combinations"""
    input = {
        "radius": radius,
        "velocity": velocity,
        "density": 3000.0,
        "strength": 1e5,
        "angle": 45.0,
        "init_altitude": 100e3,
    }
    result = planet.solve_atmospheric_entry(**input)

    # Assertions
    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0
    assert np.all(result["mass"] >= -1e-9)
    assert np.all(result["altitude"] >= -1e-6)


@mark.parametrize(
    "density,strength",
    [
        (1000.0, 1e4),  # Ice-like (weak)
        (3000.0, 1e6),  # Stony (medium)
        (5000.0, 1e7),  # Stony-iron (strong)
        (8000.0, 1e8),  # Iron (very strong)
    ],
)
def test_density_strength_combinations(planet, density, strength):
    """Test realistic density and strength combinations"""
    input = {
        "radius": 10.0,
        "velocity": 2.0e4,
        "density": density,
        "strength": strength,
        "angle": 45.0,
        "init_altitude": 100e3,
    }
    result = planet.solve_atmospheric_entry(**input)

    # Assertions
    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0
    initial_mass = (4 / 3) * np.pi * 10.0**3 * density
    assert np.isclose(result["mass"].iloc[0], initial_mass, rtol=1e-2)


@mark.parametrize(
    "angle,velocity",
    [
        (10.0, 1e4),  # Shallow, slow
        (20.0, 2e4),  # Shallow, medium
        (30.0, 3e4),  # Medium angle, fast
        (45.0, 2e4),  # 45 degrees, medium
        (45.0, 5e4),  # 45 degrees, very fast
        (60.0, 3e4),  # Steep, fast
        (75.0, 2e4),  # Very steep, medium
        (85.0, 1.5e4),  # Nearly vertical, slow
    ],
)
def test_angle_velocity_combinations(planet, angle, velocity):
    """Test realistic angle and velocity combinations"""
    input = {
        "radius": 10.0,
        "velocity": velocity,
        "density": 3000.0,
        "strength": 1e5,
        "angle": angle,
        "init_altitude": 100e3,
    }
    result = planet.solve_atmospheric_entry(**input)

    # Assertions
    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0
    assert np.isclose(result["angle"].iloc[0], angle, rtol=1e-2)
    assert np.isclose(result["velocity"].iloc[0], velocity, rtol=1e-2)
