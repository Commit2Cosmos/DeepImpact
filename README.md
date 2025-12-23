# DeepImpact

## Description

The airburst solver uses the asteroid’s physical properties and entry conditions to simulate how it behaves as it travels through the atmosphere. From these inputs, it generates a detailed trajectory profile that includes the object’s position, speed, mass, size, pancake spreading behaviour occurs, and the rate at which it loses kinetic energy. Beyond this full evolution history, the solver also identifies five essential outcomes: whether the event results in an airburst or a ground impact, the altitude where maximum energy is deposited, the maximum energy loss per kilometre, the total energy released at the burst, and the horizontal distance travelled before the burst occurs. 


## Installation

To install the module and any pre-requisites, from the base directory run
```
pip install -e .
```  

## Downloading postcode data

To download the postcode data run
```
python download_data.py
```

## dependencies

- numpy
- scipy
- sympy
- pandas
- matplotlib
- mpltools
- pytest
- sphinx

## Quick Example

### **1. Run the atmospheric entry solver**

Create a `Planet` instance and call `solve_atmospheric_entry()`.
This function computes the dynamical evolution of an object entering the atmosphere based on its initial radius, velocity, density, strength, entry angle, and other parameters. The simulation tracks:

* Changes in velocity
* Mass loss through ablation
* Remaining mass
* Evolution of the trajectory angle
* Altitude and horizontal distance
* Radius growth due to pancake-style deformation
* Spreading rate

It returns a `DataFrame` containing the complete trajectory history.

### 2.Calculate energy loss. 

Pass the perivious result about volicity and mass from the 'solve_atmospheric_enery' solution to `calculate_energy()`.

This function calculates: the rate of change of kinetic energy with altitude (dE/dz), in kt TNT/km, and adds this column to the results table.

### 3. Analyze Event Results

Use `analyze_outcome()` to analyze the results.

It identifies:
Whether the event was an airburst or a cratering;
Peak energy deposition (peak dE/dz);
Airburst altitude;
Horizontal location of the airburst;
Energy released (kt TNT).

All information is returned as a dictionary.


## Important files
To use the core functionality of the package, you can import the main components directly from the deepimpact module:

from deepimpact.solver import Planet

The Planet class provides access to the three key methods:
1. ```solve_atmospheric_entry()``` — runs the RK4 atmospheric-entry simulation
2. ```calculate_energy()``` — computes the altitude-dependent energy-loss profile
3. ```analyse_outcome()``` — identifies whether the event results in an airburst or ground impact

Example:
```
planet = Planet()
result = planet.solve_atmospheric_entry(...)
result = planet.calculate_energy(...)
outcome = planet.analyse_outcome(...)
```

To evaluate the numerical accuracy of different integration methods, you can import the RMSE comparison utilities:

from deepimpact.rk4_fwd_comp import Planet, calculate_error

The calculate_error() function computes the RMSE between the Forward Euler and RK4 solutions, and these results can be visualised through line charts and bar charts for velocity and altitude.

Plotting Energy–Altitude Profiles
To visualise how key physical quantities evolve with altitude during atmospheric entry (e.g., energy deposition, kinetic energy, trajectory), you can import the plotting utilities:

from deepimpact.plot_graphs import plot_altitude_energy

### **Numerical Choices: Solver Algorithm and Time Step (`dt`)**

The atmospheric-entry solver uses a fourth-order Runge–Kutta (RK4) integrator to evolve the governing differential equations describing velocity, mass loss, angle, altitude, radius growth, and pancaking behaviour. RK4 provides a stable, high-accuracy solution even in highly nonlinear regimes such as rapid fragmentation or strong deceleration.

The user-specified time step `dt` controls the resolution of the output. For numerical stability and accuracy, the solver internally subdivides the RK4 step:  
- The **internal RK4 step size** is constant (0.005) to ensure sufficient resolution during periods of rapid change.  
- The **output sampling step** remains equal to the user-provided `dt`.

A smaller `dt` improves accuracy in both the trajectory and energy-loss profiles, but increases computation time. Larger values of `dt` run faster but may miss sharp features in deceleration, ablation rate, or pancaking onset. For most scenarios, values of 0.005 seconds offer a good balance between precision and performance.

To prevent unnecessary computation, the solver includes an **early-termination criterion** based on the stabilisation of velocity, mass, and horizontal distance over several RK4 iterations. This allows the simulation to stop once the trajectory has converged or the object has fully ablated.


## Automated testing

To run the pytest test suite, from the base directory run
```
pytest tests/
```

Note that you should keep the tests provided, adding new ones as you develop your code. If any of these tests fail it is likely that the scoring algorithm will not work.

## Documentation

To generate the documentation (in html format)
```
python -m sphinx docs html
```

See the `docs` directory for the preliminary documentation provided that you should add to.

## Example usage

For example usage see `example.py` in the examples folder:
```
python examples/example.py
```