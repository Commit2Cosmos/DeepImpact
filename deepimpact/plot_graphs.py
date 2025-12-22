from typing import List
import matplotlib.pyplot as plt
import pandas as pd
from .utils import kinetic_energy

__all__ = ["plot_altitude_energy"]


def plot_altitude_energy(
    results: List[pd.DataFrame],
    figsize: tuple[float, float] = (18, 6),
) -> tuple[plt.Figure, tuple[plt.Axes, plt.Axes, plt.Axes]]:
    """
    Plot altitude vs. (i) dE/dz, (ii) kinetic energy, and (iii) trajectory side-by-side.

    Parameters
    ----------
    results : List[DataFrame]
        List of DataFrames from the solver. Each DataFrame should contain:
        - For energy plots: ``altitude``, ``mass``, ``velocity``, and ``dedz`` columns
        - For trajectory plot: ``altitude`` and ``distance`` columns
    figsize : tuple, optional
        Figure size passed to ``plt.subplots``. Default is (18, 6) for 3 plots.

    Returns
    -------
    fig : Figure
        The matplotlib figure object.
    axes : tuple
        Tuple of three matplotlib axes objects (ax_loss, ax_energy, ax_trajectory).
    """
    if not results:
        raise ValueError("results list cannot be empty")

    # Create 3 subplots side-by-side
    fig, axes = plt.subplots(2, 3, figsize=figsize)
    axes = axes.ravel()
    ax_loss, ax_energy, ax_trajectory, ax_angle, ax_radius, ax_mass = axes

    # Color palette for multiple results
    colors = [
        "tab:blue",
        "tab:orange",
        "tab:green",
        "tab:red",
        "tab:purple",
        "tab:brown",
    ]
    labels = ["Calculated", "True"]

    # Plot each result in the list
    for idx, result in enumerate(results):
        color = colors[idx % len(colors)]

        # Check required columns for energy plots
        required_energy = {
            "altitude",
            "mass",
            "velocity",
            "distance",
            "dedz",
            "angle",
            "radius",
        }
        missing = required_energy.difference(result.columns)
        if missing:
            raise ValueError(
                f"Input dataframe {idx} missing columns for plots: {missing}"
            )

        altitude_km = result["altitude"] / 1e3

        # Plot energy loss rate
        ax_loss.plot(result["dedz"], altitude_km, color=color, label=labels[idx])

        # Plot kinetic energy
        ax_energy.plot(
            kinetic_energy(result["mass"], result["velocity"]), altitude_km, color=color
        )

        # Plot trajectory
        distance_km = result["distance"] / 1e3
        ax_trajectory.plot(distance_km, altitude_km, color=color)

        # Plot angle
        ax_angle.plot(result["angle"], altitude_km, color=color)

        # Plot radius
        ax_radius.plot(result["radius"], altitude_km, color=color)

        # Plot mass
        ax_mass.plot(result["mass"], altitude_km, color=color)

    # Configure energy loss rate plot
    ax_loss.set_xlabel("Energy loss (kT TNT / km)")
    ax_loss.set_ylabel("Altitude (km)")
    ax_loss.set_title("Altitude vs. Energy Loss Rate")
    ax_loss.grid(True, linestyle="--", alpha=0.4)

    # Configure kinetic energy plot
    ax_energy.set_xlabel("Kinetic energy (kT)")
    ax_energy.set_ylabel("Altitude (km)")
    ax_energy.set_title("Altitude vs. Kinetic Energy")
    ax_energy.grid(True, linestyle="--", alpha=0.4)

    # Configure trajectory plot
    ax_trajectory.set_xlabel("Horizontal distance (km)")
    ax_trajectory.set_ylabel("Altitude (km)")
    ax_trajectory.set_title("Trajectory: Altitude vs. Horizontal Distance")
    ax_trajectory.grid(True, linestyle="--", alpha=0.4)

    # Configure angle plot
    ax_angle.set_xlabel("Angle (rad)")
    ax_angle.set_ylabel("Altitude (km)")
    ax_angle.set_title("Altitude vs. Angle")
    ax_angle.grid(True, linestyle="--", alpha=0.4)

    # Configure radius plot
    ax_radius.set_xlabel("Radius (m)")
    ax_radius.set_ylabel("Altitude (km)")
    ax_radius.set_title("Altitude vs. Radius")
    ax_radius.grid(True, linestyle="--", alpha=0.4)

    # Configure mass plot
    ax_mass.set_xlabel("Mass (kg)")
    ax_mass.set_ylabel("Altitude (km)")
    ax_mass.set_title("Altitude vs. Mass")
    ax_mass.grid(True, linestyle="--", alpha=0.4)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels)

    fig.tight_layout()
    plt.show()

    return fig, axes
