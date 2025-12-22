Project 1: Deep Impact - The hazard of small asteroids
=====================================================

Synopsis:
---------

Asteroids entering Earth’s atmosphere are subject to extreme drag forces
that decelerate, heat and disrupt the space rocks. The fate of an
asteroid is a complex function of its initial mass, speed, trajectory
angle and internal strength.

`Asteroids <https://en.wikipedia.org/wiki/Asteroid>`__ 10-100 m in
diameter can penetrate deep into Earth’s atmosphere and disrupt
catastrophically, generating an atmospheric disturbance
(`airburst <https://en.wikipedia.org/wiki/Air_burst>`__) that can cause
`damage on the ground <https://www.youtube.com/watch?v=tq02C_3FvFo>`__.
Such an event occurred over the city of
`Chelyabinsk <https://en.wikipedia.org/wiki/Chelyabinsk_meteor>`__ in
Russia, in 2013, releasing energy equivalent to about 520 `kilotons of
TNT <https://en.wikipedia.org/wiki/TNT_equivalent>`__ (1 kt TNT is
equivalent to :math:`4.184 \times 10^{12}` J), and injuring thousands of
people (`Popova et al.,
2013 <http://doi.org/10.1126/science.1242642>`__; `Brown et al.,
2013 <http://doi.org/10.1038/nature12741>`__). An even larger event
occurred over
`Tunguska <https://en.wikipedia.org/wiki/Tunguska_event>`__, a
relatively unpopulated area in Siberia, in 1908.

This simulator predicts the fate of asteroids entering Earth’s atmosphere,
and provides a hazard mapper for an impact over the UK.

Problem definition
------------------

Equations of motion for a rigid asteroid
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The dynamics of an asteroid in Earth’s atmosphere prior to break-up is
governed by a coupled set of ordinary differential equations:

.. math::
   :nowrap:

   \begin{aligned} 
   \frac{dv}{dt} & = \frac{-C_D\rho_a A v^2}{2 m} + g \sin \theta \\
   \frac{dm}{dt} & = \frac{-C_H\rho_a A v^3}{2 Q} \\
   \frac{d\theta}{dt} & = \frac{g\cos\theta}{v} - \frac{C_L\rho_a A v}{2 m} - \frac{v\cos\theta}{R_P + z} \\
   \frac{dz}{dt} & = -v\sin\theta \\
   \frac{dx}{dt} & = \frac{v\cos\theta}{1 + z/R_P}
   \end{aligned}

In these equations, :math:`v`, :math:`m`, and :math:`A` are the asteroid
speed (along trajectory), mass and cross-sectional area, respectively.
We will assume an initially **spherical asteroid** to convert from
inital radius to mass (and cross-sectional area). :math:`\theta` is the
meteoroid trajectory angle to the horizontal (in radians), :math:`x` is
the downrange distance of the meteoroid from its entry position,
:math:`z` is the altitude and :math:`t` is time; :math:`C_D` is the drag
coefficient, :math:`\rho_a` is the atmospheric density (a function of
altitude ), :math:`C_H` is an ablation efficiency coefficient, :math:`Q`
is the specific heat of ablation; :math:`C_L` is a lift coefficient; and
:math:`R_P` is the planetary radius. All terms use MKS units.

Asteroid break-up and deformation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A commonly used criterion for the break-up of an asteroid in the
atmosphere is when the ram pressure of the air interacting with the
asteroid :math:`\rho_a v^2` first exceeds the strength of the asteroid
:math:`Y`.

.. math:: \rho_a v^2 = Y

Should break-up occur, the asteroid deforms and spreads laterally as it
continues its passage through the atmosphere. Several models for the
spreading rate have been proposed. In this project, it will be assumed
that the fragmented asteroid’s spreading acceleration is related to its
along trajectory speed
`(Chyba et al., 1993) <https://doi.org/10.1038/361040a0>`__:

.. math::  r\frac{d^2r}{dt^2} = C_D \frac{\rho_a}{\rho_m} v^2

Where :math:`r` is the asteroid radius, :math:`\rho_m` is the asteroid
density (assumed constant) and :math:`C_D` is the drag coefficient. It
is conventional to define the cross-sectional area of the expanding
cloud of fragments as :math:`A = \pi r^2` (i.e., assuming a circular
cross-section), for use in the above equations. So, the originally
spherical asteroid spreads laterally, flattening into a "pancake".
Spreading **ceases** when the lateral expansion reaches a specified
limit :math:`r/r_0 = f_p` known as the pancake factor (:math:`r_0`
is the original radius of the asteroid), which is often assumed to
be between 4 and 7.

Airblast damage
~~~~~~~~~~~~~~~

The rapid deposition of energy in the atmosphere is analogous to an
explosion and so the environmental consequences of the airburst can be
estimated using empirical data from atmospheric explosion experiments
`(Glasstone and Dolan,
1977) <https://www.dtra.mil/Portals/61/Documents/NTPR/4-Rad_Exp_Rpts/36_The_Effects_of_Nuclear_Weapons.pdf>`__.

The main cause of damage close to the impact site is a strong (pressure)
blastwave in the air, known as the **airblast**. Empirical data suggest
that the pressure in this wave :math:`p` (in Pa) (above ambient, also
known as overpressure), as a function of explosion energy :math:`E_k`
(in kilotons of TNT equivalent), burst altitude :math:`z_b` (in m) and
horizontal range :math:`r` (in m), is given by:

.. math::
   :nowrap:

   \begin{equation*}
      p(r) = 3 \times 10^{11} \left(\frac{r^2 + z_b^2}{E_k^{2/3}}\right)^{-1.3} + 2 \times 10^{7} \left(\frac{r^2 + z_b^2}{E_k^{2/3}}\right)^{-0.57}
   \end{equation*}

For airbursts, we will take the total kinetic energy lost by the
asteroid at the burst altitude as the burst energy :math:`E_k`. For
cratering events, we will define :math:`E_k`
as the **larger** of the total kinetic energy lost by the asteroid at
the burst altitude or the residual kinetic energy of the asteroid when
it hits the ground.

The following threshold pressures can then be used to define different
degrees of damage.

+--------------+-------------------------------------+----------------+
| Damage Level | Description                         | Pressure (kPa) |
+==============+=====================================+================+
| 1            | ~10% glass windows shatter          | 1              |
+--------------+-------------------------------------+----------------+
| 2            | ~90% glass windows shatter          | 5              |
+--------------+-------------------------------------+----------------+
| 3            | Wood frame buildings collapse       | 25             |
+--------------+-------------------------------------+----------------+
| 4            | Multistory brick buildings collapse | 40             |
+--------------+-------------------------------------+----------------+

Table 1: Pressure thresholds (in kPa) for airblast damage

Additional sections
~~~~~~~~~~~~~~~~~~~

You should expand this documentation to include explanatory text for all components of your tool. 


**Damage.py:**

1.``damage_zones``
^^^^^^^^^^^^^^^^^^^

This function computes the surface zero location and the blast radii (in meters)
for a set of overpressure thresholds.

**Purpose:**

The ``damage_zones`` function turns an airburst outcome into:

- The latitude/longitude of the surface zero point
- One blast radius per requested damage threshold

It is used later by ``impact_risk`` to decide which postcodes are affected.

**Mathematical implementation:**

The ground overpressure is modelled as

.. math::

   p(r) = 3\times 10^{11} S^{-1.3} + 2\times 10^{7} S^{-0.57},

where

.. math::

   S = \frac{r^{2} + z_b^{2}}{E_k^{2/3}}.

Here :math:`r` is ground range, :math:`z_b` the burst altitude and :math:`E_k`
the burst energy in kilotons. For each target pressure :math:`p_{\mathrm{target}}`,
a bisection search finds the smallest :math:`r` with :math:`p(r) < p_{\mathrm{target}}`.


2.``impact_risk``
^^^^^^^^^^^^^^^^^^

This function performs an uncertainty analysis over many impact scenarios and
estimates both postcode risk and population affected.

**Purpose:**

Given a list of possible impacts, ``impact_risk``:

- Solves the atmospheric entry and outcome for each scenario
- Computes the blast radius at a single pressure level
- Uses the geospatial locator to find postcodes and populations inside the blast
- Aggregates these results into postcode probabilities and population statistics

**Probabilistic formulation:**

If a postcode :math:`i` is inside the blast radius :math:`n_i` times out of
:math:`N` scenarios, its estimated impact probability is

.. math::

   P_i = \frac{n_i}{N}.

The total population affected is treated as a random variable; the function
returns its sample mean and standard deviation over all scenarios.

**Airburst Solver**
 
Solver.py
1.''solve_atmospheric_entry''
 
The airburst solver integrates the physical evolution of the asteroid
according to the governing differential equations for velocity, mass,
trajectory angle, altitude, and horizontal distance. These equations
describe the interaction between the asteroid and the atmosphere before
and after fragmentation.
 
**Purpose:**
 
These equations define how the solver updates the asteroid’s physical
state during atmospheric entry. They capture the effects of drag,
blation, gravity, lift, and planetary curvature, and they include the
lateral spreading behaviour triggered when dynamic pressure exceeds the
material strength.
 
**Mathematical formulation:**
 
The solver evolves the asteroid state using the following
system of ordinary differential equations:
 
.. math::
   :nowrap:
 
   \begin{equation*}
   \frac{dv}{dt} = \frac{-C_D \rho_a A v^2}{2m} + g \sin\theta
   \end{equation*}
 
.. math::
   :nowrap:
 
   \begin{equation*}
   \frac{dm}{dt} = \frac{-C_H \rho_a A v^3}{2Q}
   \end{equation*}
 
.. math::
   :nowrap:
 
   \begin{equation*}
   \frac{d\theta}{dt} = \frac{g \cos\theta}{v} - \frac{C_L \rho_a A v}{2m} - \frac{v \cos\theta}{R_P + z}
   \end{equation*}
 
.. math::
   :nowrap:
 
   \begin{equation*}
   \frac{dz}{dt} = -v \sin\theta
   \end{equation*}
 
.. math::
   :nowrap:
 
   \begin{equation*}
   \frac{dx}{dt} = \frac{v \cos\theta}{1 + z / R_P}
   \end{equation*}
 
Fragmentation occurs when the dynamic pressure exceeds the asteroid’s strength:
 
.. math::
   :nowrap:
 
   \begin{equation*}
   \rho_a v^2 = Y
   \end{equation*}
 
After fragmentation, the asteroid spreads laterally according to:
 
.. math::
   :nowrap:
 
   \begin{equation*}
   r \frac{d^2 r}{dt^2} = C_D \frac{\rho_a}{\rho_m} v^2
   \end{equation*}
 
Spreading continues until the radius reaches the pancake limit:
 
.. math::
   :nowrap:
 
   \begin{equation*}
   r \le f_p r_0
   \end{equation*}
 
These equations form the basis of the numerical Runge–Kutta 4
integration implemented in ``solve_atmospheric_entry``.
 
**Algorithm**
 
The airburst solver integrates the governing equations using a n
umerical Runge–Kutta scheme and evaluates all physical processes
that influence the asteroid’s atmospheric entry. The algorithm
proceeds through the following steps:
 
1. **Initialisation of the state variables**  
   Using the user-provided input parameters (radius, velocity, density, s
   trength, and trajectory angle), the solver constructs the initial
   physical state of the asteroid, including mass, altitude, horizontal
    distance, radius, and lateral spreading rate.
 
2. **Atmospheric density evaluation**  
   The atmospheric density at each altitude is obtained either from an
   exponential model or from a tabulated dataset, depending on the user's
   configuration.
 
3. **Runge–Kutta time integration**  
   The coupled differential equations for velocity, mass, trajectory angle,
   altitude, horizontal distance, and radius evolution are integrated using
   a 4th-order Runge–Kutta (RK4) method.
 
4. **Fragmentation assessment**  
   At every timestep, the solver evaluates whether the dynamic
   pressure exceeds the material strength. If fragmentation occurs,
   the lateral spreading model is activated. Spreading acceleration
   continues until the radius reaches the specified pancake factor.
 
5. **State sampling**  
   To produce a manageable output size, only selected timesteps are
   recorded. Each saved state contains velocity, mass, angle, altitude,
   distance, radius, spreading rate, and time.
 
6. **Stopping conditions**  
   The integration terminates when one of the following occurs:
   - the asteroid reaches the ground,
   - the mass reaches zero (complete ablation),
   - the system enters a dynamically steady regime, where changes
   in velocity, mass, and horizontal distance fall below specified tolerances.
 
7. **Energy computation**  
   After the trajectory is produced, kinetic energy is calculated at each
   timestep, and the derivative with respect to altitude is used to obtain
   the energy-loss rate.
 
8. **Outcome determination**  
   The solver identifies the altitude at which the energy-loss rate is
   maximised and classifies the event as an airburst or a cratering impact.
 
 
2. **calculate_energy**
 
The energy calculation module computes the rate at which the asteroid
loses kinetic energy with altitude, producing the ``dedz`` profile used
to identify the burst altitude and classify the event.
 
**Purpose**
 
The purpose of ``calculate_energy`` is to derive the kinetic energy
at each recorded timestep of the trajectory, compute its rate of change
with respect to altitude, and generate a quantity (energy loss per kilometre)
that can be used to locate the peak energy deposition point.  
This value is essential for distinguishing between airburst and
cratering events and for subsequent hazard assessment.
 
**Algorithm**
 
The procedure follows these steps:
 
1. **Kinetic energy evaluation**  
   For each row of the trajectory dataset, the kinetic energy is computed
   from the instantaneous mass and velocity.
 
2. **Unit conversion**  
   The kinetic energy is converted from joules to kilotons TNT for
   consistency with hazard mapping.
 
3. **Gradient calculation**  
   Numerical gradients of the energy and altitude are computed
   with respect to time.
 
4. **Energy-loss rate**  
   The solver divides the energy gradient by the altitude gradient, yielding
   the rate of kinetic energy lost per metre of descent.
 
5. **Conversion to per-kilometre scale**  
   The energy-loss rate is multiplied by 1000 to express the final
   quantity in kilotons TNT per kilometre.
 
6. **Appending results**  
   A new column, ``dedz``, is added to the DataFrame and returned.
 
3. **analyse_outcome**
 
### **Algorithm**
 
The outcome-analysis procedure evaluates the completed atmospheric-entry solution and determines whether the event produces an airburst or a ground impact. The steps are:
 
1. **Identify the energy-loss profile**  
   Using the previously computed `dedz` values (kinetic-energy loss per kilometre of descent), the solver examines how the energy deposition varies with altitude.
 
2. **Peak detection**  
   A peak-finding algorithm is applied to the `dedz` curve. Small numerical fluctuations are filtered out by requiring a minimum prominence proportional to the global maximum.  
   A physically meaningful airburst corresponds to the most prominent peak in this energy-loss profile.
 
3. **Airburst classification**  
   If a significant peak is present, the entry is classified as an **Airburst**.  
   The solver extracts four key quantities at the peak location:  
   - maximum energy-loss rate (`burst_peak_dedz`)  
   - altitude at peak energy deposition (`burst_altitude`)  
   - down-range distance at the burst (`burst_distance`)  
   - instantaneous kinetic energy at the peak (`burst_energy`)
 
4. **Ground-impact classification**  
   If no valid peak is found, the event is designated as **Cratering** (ground impact).  
   In this case, the solver uses the following quantities at the final time step:  
   - final energy-loss rate  
   - altitude = 0  
   - horizontal impact distance  
   - the larger of:  
     • the residual kinetic energy at ground impact  
     • the total kinetic energy lost during atmospheric entry
 
5. **Outcome packaging**  
   All derived quantities are stored in a dictionary and returned, enabling downstream visualisation, damage estimation, or mapping.
 


Locator.py
~~~~~~~~~~

1.``great_circle_distance``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This function calculates the great circle distance (in meters) between pairs of points specified 
by latitude and longitude on a spherical Earth (with a radius of 6371 km).

**Purpose:**

The ``great_circle_distance`` function is repeatedly called in the algorithm for finding postcodes 
and populations within a specified zone around an impact site. Accurate distance calculation is 
essential for:

- Determining which postcodes fall within damage zones
- Identifying the nearest census grid points to an impact location

**Mathematical Implementation:**

We implement this using the **spherical law of cosines**:

.. math::
   
   \frac{r}{R_p} = \arccos\left(\sin\varphi_1 \sin\varphi_2 + \cos\varphi_1 \cos\varphi_2 \cos|\lambda_1 - \lambda_2|\right)

Where:

- :math:`r` : The great circle distance between the two points
- :math:`R_p` : Earth's radius
- :math:`\varphi_1, \varphi_2` : The latitude of point 1 and point 2
- :math:`\lambda_1, \lambda_2` : The longitude of point 1 and point 2

2. ``get_postcodes_by_radius``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We retrieves all postcodes and their corresponding latitude and longitude data from the 
``full_postcodes.csv`` file.

**Purpose:**

Given a damage center X and radii of different damage levels, the function returns postcodes 
within specific distances of the input location.

**Algorithm:**

1. Read the ``full_postcodes.csv`` file containing postcode location data
2. For each specified radius value:
   
   - Use the ``great_circle_distance`` function (mentioned in the previous section) to calculate 
     the great circle distance between each postcode and the damage center X
   - Filter out the postcodes within the blast radius by comparing calculated distances with 
     the specified radius
   - Return the list of postcodes that fall within the damage zone

3. ``get_population_by_radius``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``UK_residential_population_2011_latlon.asc`` file divides the UK into 1km × 1km grids, 
providing the center coordinates and corresponding population count for each grid.

**Purpose:**

Given a damage center X and radii of different damage levels, the function returns the population 
within specific distances of the input location.

**Algorithm:**

To locate the k nearest grid centers from damage center X, we use **progressive spatial filtering**:

- Expand search radius until enough points are found
- Compute corresponding distances using ``great_circle_distance``

For different radii, identify varying numbers of grid centers closest to X (adaptive selection based on radius size; see code for details)

For each grid, we consider four scenarios of its relationship with the damage zone to estimate 
the population:

1. The grid is inside the damage zone
2. The damage zone is inside the grid  
3. They intersect with each other
4. They are separate from each other

.. figure:: ../explain/damage_zone_grid_cases.png
   :width: 600px
   :align: center
   :alt: Four scenarios of grid-damage zone intersection
   
   Figure: Four scenarios showing the relationship between census grids (1km × 1km) and 
   circular damage zones, used for estimating affected population.

For case 1, we need to take the population of multiple grids and sum them to estimate the population of the damage zone (the number of grids to take can be dynamically selected based on the size of the radius).

For case 2, we estimate by multiplying the grid's population by the area ratio of the damage zone to the grid.

For case 3, we estimate using a linear combination of the radius, grid center, and the grid's semi-diagonal length.

For case 4, the population is 0 for the currently considered grid.



Mapping and Visualisation Tools
-------------------------------

In addition to the physical modelling and atmospheric entry solver, this project
includes a suite of interactive geospatial mapping tools designed to visualise
hazard footprints, postcode-level impact probabilities, and evacuation routes
for use in emergency planning scenarios. These visualisations transform raw
model outputs into intuitive geographic representations that support rapid
decision-making.

All maps are generated using the ``folium`` library and exported as standalone
HTML files inside the ``ui_output/`` directory. They can be viewed directly in
a web browser without requiring the Python environment.

Atmospheric Damage Map (Blast Radius Visualisation)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using results from the atmospheric entry solver and airburst damage model,
this map plots four concentric damage zones corresponding to standard
overpressure thresholds:

1 kPa  – light window damage
5 kPa  – severe window damage
25 kPa – structural damage
40 kPa – collapse of brick buildings

Each threshold is visualised using a coloured circular region to indicate the
approximate ground area affected by the airburst. This provides a clear,
interpretable link between physical model predictions and real-world
environmental impact.

To generate the damage map, run:

``python mapping_interface.py``

This script automatically loads the first scenario from
``resources/impact_parameter_list.csv`` and exports the visualisation to
``ui_output/damage_map.html``.


The output is saved as:

``ui_output/damage_map.html``

Heatmap Representation of Probabilistic Impact Footprint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The heatmap visualisation provides a smoothed, continuous representation of
postcode-level risk. Instead of discrete circles, risk values are aggregated
to form a gradient field that reveals:

* the geographical extent of the hazard footprint,
* central regions of consistently high predicted impact,
* diffuse outer areas with lower probability.

This format is well suited for communicating broad hazard patterns to
policymakers or non-technical audiences.

The heatmap is generated within the notebook
``impact_risk_visualization.ipynb`` and exported to:

``ui_output/impact_risk_heatmap.html``.


Impact Risk Map (Circle-Marker Representation)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This visualisation uses ``folium.CircleMarker`` to display postcode-level
impact probabilities. Each postcode is plotted as a circular marker whose size
and colour encode the estimated probability of experiencing damaging airblast
effects. Larger and darker markers represent higher risk.

This map highlights:

* clusters of consistently high-risk postcodes,
* gradients between high- and low-probability regions, and
* uncertainty in outer, low-probability zones.

The Circle-Marker map is generated within the notebook
``impact_risk_visualization.ipynb`` and exported to:

``ui_output/impact_risk_circles.html``

Escape Route Planner (OSRM Routing)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To demonstrate the operational use of the model in emergency response, the
escape-route module computes a path from a user’s postcode to the nearest
low-risk area. This tool:

* identifies the user’s postcode and associated risk value,
* locates the nearest postcode with probability ``p < 0.25`` (safe zone),
* queries the OSRM public API to generate a road-network route,
* overlays the route on top of the risk map for contextual awareness.

The resulting map displays the user’s location, the safe-zone marker, and the
complete OSRM route, allowing emergency planners to visualise navigable escape
paths under a hazardous impact scenario.

The resulting map is generated within the notebook
``impact_risk_visualization.ipynb`` and exported to:

``ui_output/escape_route_map.html``

Relevance for Emergency Planning
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These mapping components serve as the visual interface between the
DeepImpact model and real-world hazard interpretation. They allow emergency
management personnel to:

* assess spatial risk rapidly,
* identify heavily impacted regions,
* communicate hazard zones to stakeholders,
* and plan safe, navigable escape routes.

By integrating physical simulation outputs with real geographic data, the
mapping tools demonstrate the practical applicability of the DeepImpact
framework in public safety, contingency planning, and hazard communication.



Function API
============

.. automodule:: deepimpact
  :members:
  :imported-members:
