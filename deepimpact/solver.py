"""
This module contains the atmospheric entry solver class
for the Deep Impact project
"""

import os
import numpy as np
import pandas as pd
from .utils import joules2kilotons, kinetic_energy
from scipy.interpolate import interp1d

__all__ = ["Planet"]


class Planet:
    """
    The class called Planet is initialised with constants appropriate
    for the given target planet, including the atmospheric density profile
    and other constants
    """

    def __init__(
        self,
        atmos_func="exponential",
        atmos_filename=os.sep.join(
            (os.path.dirname(__file__), "..", "resources", "AltitudeDensityTable.csv")
        ),
        Cd=1.0,
        Ch=0.1,
        Q=1e7,
        Cl=1e-3,
        alpha=0.3,
        Rp=6371e3,
        g=9.81,
        H=8000.0,
        rho0=1.2,
        pancake_factor=7.0,
    ):
        """
        Set up the initial parameters and constants for the target planet

        Parameters
        ----------
        atmos_func : string, optional
            Function which computes atmospheric density, rho, at altitude, z.
            Default is the exponential function rho = rho0 exp(-z/H).
            Options are 'exponential', 'tabular' and 'constant'

        atmos_filename : string, optional
            Name of the filename to use with the tabular atmos_func option

        Cd : float, optional
            The drag coefficient

        Ch : float, optional
            The heat transfer coefficient

        Q : float, optional
            The heat of ablation (J/kg)

        Cl : float, optional
            Lift coefficient

        alpha : float, optional
            Dispersion coefficient

        Rp : float, optional
            Planet radius (m)

        rho0 : float, optional
            Air density at zero altitude (kg/m^3)

        g : float, optional
            Surface gravity (m/s^2)

        H : float, optional
            Atmospheric scale height (m)

        pancake_factor : float, optional
            The factor of the initial radius at which
            fragmentation will stop. Default is 7.
        """

        # Input constants
        self.Cd = Cd
        self.Ch = Ch
        self.Q = Q
        self.Cl = Cl
        self.alpha = alpha
        self.Rp = Rp
        self.g = g
        self.H = H
        self.rho0 = rho0
        self.pancake_factor = pancake_factor
        self.atmos_filename = atmos_filename

        try:
            # set function to define atmospheric density
            if atmos_func == "exponential":
                self.rhoa = lambda z: rho0 * np.exp(-z / H)
            elif atmos_func == "tabular":
                alt_rho = pd.read_csv(
                    os.sep.join(
                        (
                            os.path.dirname(__file__),
                            "..",
                            "resources",
                            "AltitudeDensityTable.csv",
                        )
                    ),
                    sep=r"\s+",
                    skiprows=1,
                    names=["altitude", "rho"],
                )
                alt_rho["log_rho"] = np.log(alt_rho["rho"])
                interp_log = interp1d(
                    alt_rho["altitude"],
                    alt_rho["log_rho"],
                    kind="linear",
                    fill_value="extrapolate",
                )

                def density(h):
                    h = float(h)
                    return float(np.exp(interp_log(h)))

                self.rhoa = lambda z: density(z)
            elif atmos_func == "constant":
                self.rhoa = lambda x: rho0
            else:
                raise NotImplementedError(
                    "atmos_func must be 'exponential', 'tabular' or 'constant'"
                )
        except NotImplementedError:
            print("atmos_func {} not implemented yet.".format(atmos_func))
            print("Falling back to constant density atmosphere for now")
            self.rhoa = lambda x: rho0

    def solve_atmospheric_entry(
        self,
        radius,
        velocity,
        density,
        strength,
        angle,
        init_altitude=100e3,
        dt=0.005,
        output_dt=0.05,
        radians=False,
    ):
        """
        Solve the system of differential equations for a given impact scenario
        using Runge-Kutta 4 (RK4) integration.

        Parameters
        ----------
        radius : float
            The radius of the asteroid in meters

        velocity : float
            The entry speed of the asteroid in meters/second

        density : float
            The density of the asteroid in kg/m^3

        strength : float
            The strength of the asteroid in Pa

        angle : float
            The initial trajectory angle of the asteroid to the horizontal
            By default, input is in degrees. If 'radians' is set to True, the
            input should be in radians

        init_altitude : float, optional
            Initial altitude in m

        dt : float, optional
            The timestep used in the integrator

        output_dt : float, optional
            The base output timestep, in s. To limit output size to "a few rows",
            the actual output timestep will be set to 100 * dt.

        radians : logical, optional
            Whether angles should be given in degrees or radians. Default=False

        Returns
        -------
        Result : DataFrame
            A pandas dataframe containing the solution to the system.
        """

        # Angles
        if not radians:
            angle_rad = np.radians(angle)
            angle_out = angle
        else:
            angle_rad = angle
            angle_out = np.degrees(angle)

        # Initial conditions
        r0 = radius
        mass = density * (4 / 3) * np.pi * r0**3

        # State vector
        y = np.array(
            [velocity, mass, angle_rad, init_altitude, 0.0, r0, 0.0], dtype=float
        )

        t = 0.0

        # RK4 internal timestep
        next_output_time = output_dt

        # Results store
        results = {
            "velocity": [velocity],
            "mass": [mass],
            "angle": [angle_out],
            "altitude": [init_altitude],
            "distance": [0.0],
            "radius": [radius],
            "spreading_rate": [0.0],
            "time": [t],
        }

        # Early-stopping settings
        # steady_steps represents number of RK4 steps to check for steady state
        steady_time_duration = 5  # seconds
        max_steady_steps_limit = 5000

        steady_steps = min(max_steady_steps_limit, int(steady_time_duration / dt))
        rel_tol_trajectory = 5e-5
        rel_tol_dynamics = 5e-5  # 0.1% velocity and mass
        eps = 1e-12

        vel_hist = np.empty(steady_steps + 1)
        mass_hist = np.empty(steady_steps + 1)
        dist_hist = np.empty(steady_steps + 1)

        count = 0
        idx = 0

        # Constants
        Cd, Ch, Q, Cl, g = self.Cd, self.Ch, self.Q, self.Cl, self.g
        Rp = self.Rp
        rhoa_func = self.rhoa
        rho_m = density
        met_strength = strength
        Pf = self.pancake_factor

        # Derivatives
        def compute_derivatives(state):
            v_curr, m_curr, theta_curr, z_curr, x_curr, r_curr, u_curr = state

            """
            Compute the derivatives of the state vector for atmospheric entry dynamics.
        
            This function calculates the time derivatives of velocity, mass, angle, altitude,
            distance, radius, and spreading rate for an asteroid undergoing atmospheric entry.
            It includes effects of drag, ablation, lift, gravity, planetary curvature, and
            pancaking when ram pressure exceeds the asteroid's strength.
            
            Parameters
            ----------
            time : float
                The current time in seconds (not used in calculation but required for 
                compatibility with RK4 integration)
            
            state : array
                A 7-element state vector containing:
                [0] v_curr : velocity (m/s)
                [1] m_curr : mass (kg)
                [2] theta_curr : trajectory angle (radians)
                [3] z_curr : altitude (m)
                [4] x_curr : horizontal distance traveled (m)
                [5] r_curr : current radius (m)
                [6] u_curr : spreading rate (m/s)
            
            Returns
            -------
            derivatives : ndarray
                A 7-element array containing the time derivatives:
                [0] dv/dt : rate of change of velocity (m/s²)
                [1] dm/dt : rate of change of mass (kg/s)
                [2] dtheta/dt : rate of change of angle (rad/s)
                [3] dz/dt : rate of change of altitude (m/s)
                [4] dx/dt : rate of change of horizontal distance (m/s)
                [5] dr/dt : rate of change of radius (m/s)
                [6] du/dt : rate of change of spreading rate (m/s²)
            
            Notes
            -----
            The function returns zeros if the asteroid is below ground level (z <= 0)
            or has no remaining mass (m <= 0). Pancaking occurs when ram pressure
            exceeds the meteoroid strength and the current radius is less than or
            equal to Pf times the initial radius.
            """

            # Stop if below ground or no mass
            if z_curr <= 0.0 or m_curr <= 0.0:
                return np.zeros(7, dtype=float)

            # Atmosphere
            rho_a = rhoa_func(z_curr)
            ram_p = rho_a * v_curr**2

            # Cross-sectional area
            area = np.pi * r_curr**2

            # dv/dt
            dv_dt = -(Cd * rho_a * area * v_curr**2) / (2.0 * m_curr) + g * np.sin(
                theta_curr
            )

            # dm/dt (ablation)
            dm_dt = -(Ch * rho_a * area * v_curr**3) / (2.0 * Q)

            # dtheta/dt (lift + gravity + curvature)
            dtheta_dt = (
                (g * np.cos(theta_curr)) / v_curr
                - (Cl * rho_a * area * v_curr) / (2.0 * m_curr)
                - (v_curr * np.cos(theta_curr)) / (Rp + z_curr)
            )

            # dz/dt
            dz_dt = -v_curr * np.sin(theta_curr)

            # dx/dt (curved planet)
            dx_dt = v_curr * np.cos(theta_curr) / (1.0 + z_curr / Rp)

            # Pancaking when ram pressure exceeds strength and radius < Pf * r0
            is_pancaking = (ram_p >= met_strength) and (r_curr <= Pf * r0)

            if is_pancaking:
                # Pancaking equations
                du_dt = (Cd * rho_a * v_curr**2) / (rho_m * r_curr)
                dr_dt = u_curr
            else:
                # No pancaking
                du_dt = 0.0
                dr_dt = 0.0

            return np.array(
                [dv_dt, dm_dt, dtheta_dt, dz_dt, dx_dt, dr_dt, du_dt], dtype=float
            )

        # RK4 loop
        iteration = 0
        # maximum iterations cap depending on the internal timestep
        max_iter = (1000 / dt) + 5000

        while (y[3] > 0.0) and (y[1] > 0.0) and (iteration < max_iter):
            iteration += 1

            k1 = compute_derivatives(y)
            k2 = compute_derivatives(y + 0.5 * dt * k1)
            k3 = compute_derivatives(y + 0.5 * dt * k2)
            k4 = compute_derivatives(y + dt * k3)

            y_new = y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
            t_new = t + dt

            # Boundary handling (ground or full ablation)
            is_ground_hit = y_new[3] <= 0.0
            is_ablated = y_new[1] <= 0.0

            if is_ground_hit or is_ablated:
                if is_ground_hit and k1[3] < 0.0:
                    dt_remain = -y[3] / k1[3]
                    if 1e-10 < dt_remain < dt:
                        y = y + k1 * dt_remain
                        t = t + dt_remain
                    else:
                        y = y_new
                        t = t_new
                else:
                    y = y_new
                    t = t_new

                # Clamp final state
                y[3] = max(y[3], 0.0)
                y[1] = max(y[1], 0.0)

                # Store final point
                results["velocity"].append(y[0])
                results["mass"].append(y[1])
                results["angle"].append(y[2] if radians else np.degrees(y[2]))
                results["altitude"].append(y[3])
                results["distance"].append(y[4])
                results["radius"].append(y[5])
                results["spreading_rate"].append(y[6])
                results["time"].append(t)
                break

            y = y_new
            t = t_new

            # Track history at every RK4 step for early stopping
            vel_hist[idx] = y[0]
            mass_hist[idx] = y[1]
            dist_hist[idx] = y[4]

            idx = (idx + 1) % (steady_steps + 1)
            count = min(count + 1, (steady_steps + 1))

            # Early-stopping condition based on last n RK4 steps
            if count == steady_steps + 1:
                recent_v = np.concatenate((vel_hist[idx:], vel_hist[:idx]))
                recent_m = np.concatenate((mass_hist[idx:], mass_hist[:idx]))
                recent_x = np.concatenate((dist_hist[idx:], dist_hist[:idx]))

                # find absolute changes
                dv = np.abs(np.diff(recent_v))
                dm = np.abs(np.diff(recent_m))
                dx = np.abs(np.diff(recent_x))

                # relative changes
                rel_dv = dv / np.maximum(np.abs(recent_v[1:]), eps)
                rel_dm = dm / np.maximum(np.abs(recent_m[1:]), eps)
                rel_dx = dx / np.maximum(np.abs(recent_x[1:]), eps)

                # if trajectory is not changing by a relative amount in the last n steps and mass and velocity do not change then:
                steady = (
                    np.all(rel_dx < rel_tol_trajectory)
                    and np.all(rel_dv < rel_tol_dynamics)
                    and np.all(rel_dm < rel_tol_dynamics)
                )

                if steady:
                    return pd.DataFrame(results)

            # Output at output times
            if t >= next_output_time:
                # Store current state
                results["velocity"].append(y[0])
                results["mass"].append(y[1])
                results["angle"].append(y[2] if radians else np.degrees(y[2]))
                results["altitude"].append(y[3])
                results["distance"].append(y[4])
                results["radius"].append(y[5])
                results["spreading_rate"].append(y[6])
                results["time"].append(t)

                next_output_time += output_dt

        return pd.DataFrame(results)

    def calculate_energy(self, result: pd.DataFrame):
        """
        Function to calculate the kinetic energy lost per unit altitude in
        kilotons TNT per km, for a given solution.

        Parameters
        ----------
        result : DataFrame
            A pandas dataframe with columns for the velocity, mass, angle,
            altitude, horizontal distance and radius as a function of time

        Returns : DataFrame
            Returns the dataframe with additional column ``dedz``
            (kinetic energy lost per unit altitude)

        """

        result = result.copy()

        # calculate kinetic energy
        ke = kinetic_energy(result["mass"], result["velocity"])
        ke = joules2kilotons(ke)

        dKE_dt = np.gradient(ke, edge_order=2)
        dz_dt = np.gradient(result["altitude"], edge_order=2)
        result["dedz"] = (dKE_dt / dz_dt) * 1000

        return result

    def analyse_outcome(self, result: pd.DataFrame):
        """
        Inspect a pre-found solution to calculate the impact and airburst stats

        Parameters
        ----------
        result : DataFrame
            pandas dataframe with velocity, mass, angle, altitude, horizontal
            distance, radius and dedz as a function of time

        Returns
        -------
        outcome : Dict
            dictionary with details of the impact event, which should contain
            the key:
                ``outcome`` (which should contain one of the
                following strings: ``Airburst`` or ``Cratering``),
            as well as the following 4 keys:
                ``burst_peak_dedz``, ``burst_altitude``,
                ``burst_distance``, ``burst_energy``
        """

        outcome = {}

        max_dedz_idx = np.argmax(result["dedz"])

        outcome["burst_altitude"] = result["altitude"].iloc[max_dedz_idx]

        if outcome["burst_altitude"] > 0:
            outcome["outcome"] = "Airburst"
            outcome["burst_peak_dedz"] = result["dedz"].iloc[max_dedz_idx]
            outcome["burst_distance"] = result["distance"].iloc[max_dedz_idx]
            outcome["burst_energy"] = joules2kilotons(
                kinetic_energy(
                    result["mass"].iloc[max_dedz_idx],
                    result["velocity"].iloc[max_dedz_idx],
                )
            )
        else:
            outcome["outcome"] = "Cratering"
            outcome["burst_peak_dedz"] = result["dedz"].iloc[-1]
            outcome["burst_altitude"] = 0
            outcome["burst_distance"] = result["distance"].iloc[-1]

            energy_at_impact = joules2kilotons(
                kinetic_energy(result["mass"].iloc[-1], result["velocity"].iloc[-1])
            )
            energy_init = joules2kilotons(
                kinetic_energy(result["mass"].iloc[0], result["velocity"].iloc[0])
            )

            # Burst energy is the larger of:
            # - residual kinetic energy at impact
            # - total kinetic energy lost during entry
            outcome["burst_energy"] = max(
                energy_at_impact, energy_init - energy_at_impact
            )

        return outcome
