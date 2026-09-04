import deepimpact
from deepimpact.plot_graphs import plot_altitude_energy
import time


# Initialise the Planet class
earth = deepimpact.Planet()


# Solve the atmospheric entry problem for a given set of input parameters
start = time.perf_counter()
result = earth.solve_atmospheric_entry(
    radius=35, angle=45, strength=1e7, density=3000, velocity=19e3
)
end = time.perf_counter()

plot_altitude_energy([result])

print(f"{end - start:.4f} seconds")
