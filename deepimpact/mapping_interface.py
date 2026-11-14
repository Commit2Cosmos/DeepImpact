"""
Automatically loads an impact scenario from impact_parameter_list.csv
Runs atmospheric entry solver + damage_zones
Generates 4-level blast map (yellow/orange/red/purple)
"""

import folium
import pandas as pd
import deepimpact
import webbrowser
import os
from deepimpact.mapping import plot_circle
from deepimpact.damage import damage_zones


def generate_damage_map_from_csv(csv_path):
    """
    Load the first scenario from CSV and generate damage map.

    AI-GENERATED CODE NOTICE:
    Parts of this function (including structure, map-generation logic,
    and file-saving patterns) were developed with assistance from
    ChatGPT (OpenAI). Conversation summary is included in the
    project references as required.
    """

    df = pd.read_csv(csv_path)

    # Use first scenario (row 0)
    row = df.iloc[0]

    radius = row["radius"]
    angle = row["angle"]
    strength = row["strength"]
    density = row["density"]
    velocity = row["velocity"]
    lat = row["entry latitude"]
    lon = row["entry longitude"]
    bearing = row["bearing"]

    print("\n=== Loaded Scenario from CSV ===")
    print(df.iloc[0])

    # Run solver (use Planet from deepimpact)
    planet = deepimpact.Planet()

    result = planet.solve_atmospheric_entry(
        radius=radius,
        angle=angle,
        strength=strength,
        density=density,
        velocity=velocity,
        init_altitude=100e3,
    )

    result = planet.calculate_energy(result)
    outcome = planet.analyse_outcome(result)

    print("\n=== Solver Outcome ===")
    print(outcome)

    # Damage radii
    pressures = [1e3, 5e3, 25e3, 40e3]

    blat, blon, radii = damage_zones(outcome, lat, lon, bearing, pressures)

    print("\nSurface Zero:")
    print("  lat:", blat)
    print("  lon:", blon)

    print("\nDamage Radii (meters):")
    for i, r in enumerate(radii, 1):
        print(f"  Level {i}: {r:.0f} m")

    # Draw map
    fmap = folium.Map(location=[blat, blon], zoom_start=7)

    colors = ["yellow", "orange", "red", "purple"]

    for r, c in zip(radii, colors):
        plot_circle(blat, blon, r, fmap=fmap, color=c, fill_opacity=0.15)

    return fmap


if __name__ == "__main__":
    output_dir = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "ui_output")
    )
    os.makedirs(output_dir, exist_ok=True)

    # generate map
    fmap = generate_damage_map_from_csv("../resources/impact_parameter_list.csv")
    output_file = os.path.join(output_dir, "damage_map.html")
    fmap.save(output_file)
    print(f"\nMap saved as: {output_file}")

    webbrowser.open("file://" + os.path.realpath(output_file))
