"""Module to calculate the damage and impact risk for given scenarios"""

import pandas as pd
import os
import math
from collections import Counter
import numpy as np
import deepimpact

__all__ = ["damage_zones", "impact_risk"]

EARTH_RADIUS_M = 6_371_000.0  # meters


def damage_zones(outcome, lat, lon, bearing, pressures):
    """
    Calculate the latitude and longitude of the surface zero location and the
    list of airblast damage radii (m) for a given impact scenario.

    Parameters
    ----------

    outcome: Dict
        the outcome dictionary from an impact scenario
    lat: float
        latitude of the meteoroid entry point (degrees)
    lon: float
        longitude of the meteoroid entry point (degrees)
    bearing: float
        Bearing (azimuth) relative to north of meteoroid trajectory (degrees)
    pressures: float, arraylike
        List of threshold pressures to define airblast damage levels

    Returns
    -------

    blat: float
        latitude of the surface zero point (degrees)
    blon: float
        longitude of the surface zero point (degrees)
    damrad: arraylike, float
        List of distances specifying the blast radii
        for the input damage levels

    Examples
    --------

    >>> import deepimpact
    >>> outcome = {'burst_altitude': 8e3, 'burst_energy': 7e3,
                   'burst_distance': 90e3, 'burst_peak_dedz': 1e3,
                   'outcome': 'Airburst'}
    >>> deepimpact.damage_zones(outcome, 52.79, -2.95, 135,
                                pressures=[1e3, 5e3, 25e3, 40e3])
    """

    # Replace this code with your own. For demonstration we just
    # return lat, lon and a radius of 5000 m for each pressure

    # Extract required values from outcome
    try:
        burst_altitude = float(outcome["burst_altitude"])
        burst_energy = float(outcome["burst_energy"])
        burst_distance = float(outcome["burst_distance"])
    except KeyError as exc:
        raise KeyError(
            "outcome must contain 'burst_altitude', 'burst_energy' "
            "and 'burst_distance'"
        ) from exc

    # Surface zero location (point on the ground under the burst)
    blat, blon = destination_point(lat, lon, bearing, burst_distance)

    # Normalise pressures to a list and compute one radius per threshold
    if isinstance(pressures, (int, float)):
        pressure_list = [float(pressures)]
    else:
        pressure_list = [float(p) for p in pressures]

    damrad = [
        pressure_radius(p_target, burst_altitude, burst_energy)
        for p_target in pressure_list
    ]

    return blat, blon, damrad


def destination_point(lat_deg, lon_deg, bearing_deg, distance_m):
    """
    Great-circle destination given start (lat, lon), bearing and distance.
    All angles should be in degree; distance in metres.
    """

    if lat_deg < -90.0 or lat_deg > 90.0:
        raise ValueError("Latitude must be in [-90, 90] degrees.")
    if lon_deg < -180.0 or lon_deg > 180.0:
        raise ValueError("Longitude must be in [-180, 180] degrees.")
    if abs(lat_deg) > 90:
        raise ValueError("Latitude out of range [-90, 90] degrees.")

    phi1 = math.radians(lat_deg)
    lam1 = math.radians(lon_deg)
    theta = math.radians(bearing_deg)
    delta = distance_m / EARTH_RADIUS_M

    sin_phi2 = math.sin(phi1) * math.cos(delta) + math.cos(phi1) * math.sin(
        delta
    ) * math.cos(theta)
    phi2 = math.asin(sin_phi2)

    # numerical safety
    sin_phi2 = min(1.0, max(-1.0, sin_phi2))
    phi2 = math.asin(sin_phi2)

    num = math.sin(theta) * math.sin(delta) * math.cos(phi1)
    den = math.cos(delta) - math.sin(phi1) * math.sin(phi2)
    lam2 = lam1 + math.atan2(num, den)

    # Wrap longitude to [-pi, pi]
    lam2 = (lam2 + math.pi) % (2.0 * math.pi) - math.pi

    return math.degrees(phi2), math.degrees(lam2)


def overpressure(r, z_b, E_k):
    """
    Airblast overpressure p(r) in Pa for ground range r (m),
    burst altitude z_b (m) and burst energy E_k (kilotons TNT).
    """
    # No energy -> no blast
    E_k = float(E_k)
    if E_k <= 0.0:
        return 0.0

    # Dimensionless scaled distance
    scaled = (r * r + z_b * z_b) / (E_k ** (2.0 / 3.0))

    # At the burst point (scaled == 0) treat overpressure as effectively infinite
    if scaled <= 0.0:
        return float("inf")

    term1 = 3.0e11 * (scaled**-1.3)
    term2 = 2.0e7 * (scaled**-0.57)
    return term1 + term2


def pressure_radius(p_target, z_b, E_k, tol=1.0, r_init=1_000.0):
    """
    Radius at which overpressure first falls below p_target.

    Returns 0.0 if p_target is never reached on the ground
    or if the burst energy is zero/non-physical.
    """
    # If there is no energy, there can be no blast radius
    E_k = float(E_k)
    if E_k <= 0.0:
        return 0.0

    # If even at r = 0 the pressure is already below the target, no damage zone
    if overpressure(0.0, z_b, E_k) < p_target:
        return 0.0

    # Find an upper bound where pressure < target
    r_low = 0.0
    r_high = r_init
    max_range = 10_000_000.0  # safety cap
    while overpressure(r_high, z_b, E_k) > p_target and r_high < max_range:
        r_high *= 2.0

    # Bisection between r_low and r_high (catogorize pressure)
    while (r_high - r_low) > tol:
        r_mid = 0.5 * (r_low + r_high)
        if overpressure(r_mid, z_b, E_k) >= p_target:
            r_low = r_mid
        else:
            r_high = r_mid

    return 0.5 * (r_low + r_high)


def impact_risk(
    planet,
    impact_file=os.sep.join(
        (os.path.dirname(__file__), "..", "resources", "impact_parameter_list.csv")
    ),
    pressure=30.0e3,
    nsamples=None,
):
    """
    Perform an uncertainty analysis to calculate the probability for
    each affected UK postcode and the total population affected.

    Parameters
    ----------
    planet: deepimpact.Planet instance
        The Planet instance from which to solve the atmospheric entry

    impact_file: str
        Filename of a .csv file containing the impact parameter list
        with columns for 'radius', 'angle', 'velocity', 'strength',
        'density', 'entry latitude', 'entry longitude', 'bearing'

    pressure: float
        A single pressure at which to calculate the damage zone for each impact

    nsamples: int or None
        The number of iterations to perform in the uncertainty analysis.
        If None, the full set of impact parameters provided in impact_file
        is used.

    Returns
    -------
    probability: DataFrame
        A pandas DataFrame with columns for postcode and the
        probability the postcode was inside the blast radius.
    population: dict
        A dictionary containing the mean and standard deviation of the
        population affected by the impact, with keys 'mean' and 'stdev'.
        Values are floats.
    """

    # Check that pressure looks sensible
    if not isinstance(pressure, (int, float)):
        raise TypeError("pressure must be a single numeric value in Pa.")
    pressure = float(pressure)

    # Read all scenarios from the CSV
    data = pd.read_csv(impact_file)  # "resources/impact_parameter_list.csv"

    # Decide how many scenarios to use
    if nsamples is None:
        n_use = len(data)
    else:
        if nsamples <= 0:
            raise ValueError("nsamples must be a positive integer or None.")
        n_use = min(int(nsamples), len(data))
    data = data.iloc[:n_use]

    if n_use == 0:
        # No scenarios, return empty results
        probability_df = pd.DataFrame(columns=["postcode", "probability"])
        population_stats = {"mean": 0.0, "stdev": 0.0}
        return probability_df, population_stats

    # Prepare the postcode / population locator
    locator = deepimpact.GeospatialLocator()

    all_postcodes = []  # flattened list of all postcodes hit over all scenarios
    populations_per_scenario = (
        []
    )  # total population affected per scenario (including zeros)

    for i in range(n_use):
        # Extract parameters
        radius = data.loc[i, "radius"]
        angle = data.loc[i, "angle"]
        velocity = data.loc[i, "velocity"]
        strength = data.loc[i, "strength"]
        density = data.loc[i, "density"]
        entry_lat = data.loc[i, "entry latitude"]
        entry_lon = data.loc[i, "entry longitude"]
        bearing = data.loc[i, "bearing"]

        # 1. Solve atmospheric entry and outcome
        result = planet.solve_atmospheric_entry(
            radius=radius,
            angle=angle,
            strength=strength,
            density=density,
            velocity=velocity,
        )
        result = planet.calculate_energy(result)
        outcome = planet.analyse_outcome(result)

        # 2. Compute damage zone at the specified pressure
        blat, blon, radii = damage_zones(
            outcome,
            lat=entry_lat,
            lon=entry_lon,
            bearing=bearing,
            pressures=[pressure],
        )

        if len(radii) == 0 or radii[0] <= 0.0:
            # No blast radius at this pressure -> population = 0, no postcodes hit
            populations_per_scenario.append(0.0)
            continue

        damage_radius = radii[0]

        # 3. Use locator to get postcodes and population inside this radius
        postcodes_list = locator.get_postcodes_by_radius((blat, blon), [damage_radius])
        population_list = locator.get_population_by_radius(
            (blat, blon), [damage_radius]
        )

        if len(postcodes_list) == 0 or len(population_list) == 0:
            # Locator returned nothing -> treat as zero population
            populations_per_scenario.append(0.0)
            continue

        postcodes_in_radius = postcodes_list[0]
        population_in_radius = population_list[0]

        # Aggregate results
        all_postcodes.extend(postcodes_in_radius)
        populations_per_scenario.append(float(population_in_radius))

    # If no postcodes were ever hit, return empty probability but valid stats
    if len(all_postcodes) == 0:
        probability_df = pd.DataFrame(columns=["postcode", "probability"])
        mean_pop = (
            float(np.mean(populations_per_scenario))
            if populations_per_scenario
            else 0.0
        )
        std_pop = (
            float(np.std(populations_per_scenario)) if populations_per_scenario else 0.0
        )
        population_stats = {"mean": mean_pop, "stdev": std_pop}
        return probability_df, population_stats

    # Count how many times each postcode was inside the blast radius
    counts = Counter(all_postcodes)
    postcodes_sorted = sorted(counts.keys())
    probabilities = [counts[pc] / float(n_use) for pc in postcodes_sorted]

    probability_df = pd.DataFrame(
        {"postcode": postcodes_sorted, "probability": probabilities}
    )

    # Population statistics over all scenarios (including zeros)
    mean_pop = (
        float(np.mean(populations_per_scenario)) if populations_per_scenario else 0.0
    )
    std_pop = (
        float(np.std(populations_per_scenario)) if populations_per_scenario else 0.0
    )
    population_stats = {"mean": mean_pop, "stdev": std_pop}

    return probability_df, population_stats
