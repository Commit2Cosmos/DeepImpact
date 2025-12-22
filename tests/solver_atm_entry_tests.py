"""
Test file for the atmospheric entry solver
"""

import deepimpact

if __name__ == "__main__":
    planet = deepimpact.Planet(atmos_func="exponential", rho0=1.2, H=8000.0, Cd=1.0)

    result = planet.solve_atmospheric_entry(
        radius=35,
        angle=45,
        strength=1e7,
        density=3000,
        velocity=19e3,
        init_altitude=100e3,
        dt=0.1,
        radians=False,
    )

    print("Result DataFrame:")
    print(result.head())
    print(f"\nShape: {result.shape}")
    print(f"\nLast few rows:")
    print(result.tail())
    print(f"\nSimulation time: {result['time'].iloc[-1]:.2f} seconds")
    print(f"Final altitude: {result['altitude'].iloc[-1]:.2f} m")
    print(f"Final velocity: {result['velocity'].iloc[-1]:.2f} m/s")
