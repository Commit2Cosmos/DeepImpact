import numpy as np
import pandas as pd
import pytest

from deepimpact.damage import (
    damage_zones,
    destination_point,
    pressure_radius,
    impact_risk,
    overpressure,
)


# =========================
# Fixtures & Test Data
# =========================


@pytest.fixture
def outcome_airburst():
    """A representative airburst outcome used in several tests."""
    return {
        "burst_altitude": 8_000.0,  # m
        "burst_energy": 7e3,  # kilotons TNT
        "burst_distance": 90_000.0,  # m
    }


@pytest.fixture
def standard_pressures():
    """Standard four damage thresholds: 1, 5, 25, 40 kPa."""
    return [1e3, 5e3, 25e3, 40e3]


@pytest.fixture
def entry_location():
    """Canonical entry location and bearing."""
    return dict(lat=52.79, lon=-2.95, bearing=135.0)


class DummyPlanet:
    """Stub Planet that records calls and returns a fixed outcome."""

    def __init__(self):
        self.called_with = []

    def solve_atmospheric_entry(
        self, radius, velocity, density, strength, angle, **kwargs
    ):
        self.called_with.append(
            dict(
                radius=radius,
                velocity=velocity,
                density=density,
                strength=strength,
                angle=angle,
            )
        )
        # minimal dataframe with two steps so calculate_energy won't crash
        return pd.DataFrame(
            {
                "velocity": [velocity, velocity * 0.5],
                "mass": [1e5, 1e5],
                "angle": [angle, angle],
                "altitude": [100e3, 0.0],
                "distance": [0.0, 50e3],
                "radius": [radius, radius],
                "spreading_rate": [0.0, 0.0],
                "time": [0.0, 10.0],
            }
        )

    def calculate_energy(self, result):
        # Let tests ignore real energy calculation
        return result

    def analyse_outcome(self, result):
        return {
            "burst_altitude": 8_000.0,
            "burst_energy": 7e3,
            "burst_distance": 90_000.0,
            "burst_peak_dedz": 1e3,
            "outcome": "Airburst",
        }


# =========================
# damage_zones tests
# =========================


def test_damage_zones_returns_array_like_with_four_radii(
    outcome_airburst, standard_pressures, entry_location
):
    """
    damage_zones should return an array-like of 4 radii
    when given 4 pressure thresholds.
    """
    blat, blon, damrad = damage_zones(
        outcome_airburst,
        lat=entry_location["lat"],
        lon=entry_location["lon"],
        bearing=entry_location["bearing"],
        pressures=standard_pressures,
    )

    # Types for lat/lon
    assert isinstance(blat, float)
    assert isinstance(blon, float)

    # damrad should be array-like
    assert hasattr(damrad, "__len__")

    # It should contain exactly four entries (one per pressure)
    assert len(damrad) == 4
    assert len(damrad) == len(standard_pressures)

    # The container type should be something sensible
    assert isinstance(damrad, (list, tuple, np.ndarray))

    # Each entry should be a numeric scalar (float or int)
    for r in damrad:
        assert np.isscalar(r)


def test_damage_zones_monotonic_radii(
    outcome_airburst, standard_pressures, entry_location
):
    """Higher pressure thresholds should give smaller or equal radii."""
    _, _, damrad = damage_zones(
        outcome_airburst,
        lat=entry_location["lat"],
        lon=entry_location["lon"],
        bearing=entry_location["bearing"],
        pressures=standard_pressures,
    )

    # As pressure increases, radius should not increase
    for i in range(1, len(damrad)):
        assert damrad[i] <= damrad[i - 1] + 1e-6


def test_damage_zones_missing_outcome_keys_raises(standard_pressures, entry_location):
    """damage_zones should raise KeyError if outcome is incomplete."""
    outcome_incomplete = {
        "burst_altitude": 8_000.0,
        # "burst_energy" missing
        "burst_distance": 90_000.0,
    }
    with pytest.raises(KeyError):
        damage_zones(
            outcome_incomplete,
            entry_location["lat"],
            entry_location["lon"],
            entry_location["bearing"],
            standard_pressures[:2],
        )


def test_damage_zones_raises_keyerror_for_missing_outcome_keys():
    from deepimpact.damage import damage_zones

    bad_outcome = {"burst_altitude": 1e3}  # missing burst_energy, burst_distance
    with pytest.raises(KeyError):
        damage_zones(bad_outcome, 52.0, 0.0, 90.0, pressures=[1e3])


def test_damage_zones_scalar_pressure():
    from deepimpact.damage import damage_zones

    outcome = {
        "burst_altitude": 8e3,
        "burst_energy": 7e3,
        "burst_distance": 90e3,
    }
    blat, blon, damrad = damage_zones(outcome, 52.0, 0.0, 90.0, pressures=1e3)
    assert isinstance(damrad, list)
    assert len(damrad) == 1


# =========================
# pressure_radius tests
# =========================


def test_pressure_radius_too_high_returns_zero():
    """If the target pressure is unrealistically high, radius should be zero."""
    z_b = 8_000.0
    E_k = 7e3

    r = pressure_radius(1e9, z_b, E_k)  # 1 GPa, way too high
    assert r == 0.0


def test_destination_point_invalid_latitude():
    with pytest.raises(ValueError):
        destination_point(100.0, 0.0, 0.0, 1000.0)


def test_destination_point_invalid_longitude():
    with pytest.raises(ValueError):
        destination_point(0.0, 200.0, 0.0, 1000.0)


def test_overpressure_zero_energy():
    assert overpressure(1000.0, 1000.0, 0.0) == 0.0


def test_pressure_radius_zero_energy():
    assert pressure_radius(1e5, 1000.0, 0.0) == 0.0


def test_pressure_radius_no_damage_zone():
    # Overpressure at r=0 already below p_target -> radius 0
    r = pressure_radius(1e12, z_b=1e5, E_k=1.0)  # absurdly high pressure
    assert r == 0.0


# =========================
# destination_point tests
# =========================


def test_destination_point_zero_distance_same_location():
    """Zero distance should return the original coordinates."""
    lat, lon = 52.0, -3.0
    blat, blon = destination_point(lat, lon, 45.0, 0.0)

    assert pytest.approx(lat, rel=0, abs=1e-10) == blat
    assert pytest.approx(lon, rel=0, abs=1e-10) == blon


# =========================
# impact_risk tests
# =========================


def test_impact_risk_reads_csv_columns(tmp_path):
    """impact_risk should read the CSV and pass parameters into Planet.solve_atmospheric_entry."""
    # 1) make a tiny custom impact file
    csv_path = tmp_path / "impact_small.csv"
    df = pd.DataFrame(
        {
            "radius": [50.0],
            "angle": [30.0],
            "velocity": [20e3],
            "strength": [1e6],
            "density": [3_000.0],
            "entry latitude": [52.0],
            "entry longitude": [-1.5],
            "bearing": [120.0],
        }
    )
    df.to_csv(csv_path, index=False)

    planet = DummyPlanet()

    # 2) call impact_risk with our custom file
    probability, population = impact_risk(
        planet,
        impact_file=str(csv_path),
        pressure=30e3,
        nsamples=1,
    )

    # We don't care about probability/population values here, only that
    # the Planet was called correctly.
    assert len(planet.called_with) == 1
    params = planet.called_with[0]

    # Just check that these keys exist and are finite numbers
    for key in ["radius", "angle", "velocity", "strength", "density"]:
        assert key in params
        assert np.isfinite(params[key])
