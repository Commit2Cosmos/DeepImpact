"""Module dealing with postcode information."""

import numpy as np
import pandas as pd
import os

__all__ = ["GeospatialLocator", "great_circle_distance"]


def great_circle_distance(latlon1, latlon2):
    """
    Calculate the great circle distance (in metres) between pairs of
    points specified as latitude and longitude on a spherical Earth
    (with radius 6371 km).

    Parameters
    ----------

    latlon1: arraylike
        latitudes and longitudes of first point (as [n, 2] array for n points)
    latlon2: arraylike
        latitudes and longitudes of second point (as [m, 2] array for m points)

    Returns
    -------

    numpy.ndarray
        Distance in metres between each pair of points (as an n x m array)

    Examples
    --------

    >>> import numpy
    >>> fmt = lambda x: numpy.format_float_scientific(x, precision=3)}
    >>> with numpy.printoptions(formatter={'all', fmt}):
        print(great_circle_distance([[54.0, 0.0], [55, 0.0]], [55, 1.0]))
    [1.286e+05 6.378e+04]
    """
    # Ensure latlon1 and latlon2 are at least 2-dimensional
    latlon1, latlon2 = np.atleast_2d(latlon1), np.atleast_2d(latlon2)

    # Validate latitude and longitude values
    for latlon in [latlon1, latlon2]:
        if np.any(latlon[:, 0] < -90) or np.any(latlon[:, 0] > 90):
            raise ValueError("Latitude must be in the range -90 to 90")
        if np.any(latlon[:, 1] < -180) or np.any(latlon[:, 1] > 180):
            raise ValueError("Longitude must be in the range -180 to 180")

    # Convert to radians and extract coordinates
    lat1, lon1 = np.radians(latlon1).T[:, :, None]
    lat2, lon2 = np.radians(latlon2).T[:, None, :]

    # Spherical law of cosines
    cos_d = np.sin(lat1) * np.sin(lat2) + np.cos(lat1) * np.cos(lat2) * np.cos(
        lon2 - lon1
    )
    # Clipping to avoid precision issues
    cos_d = np.clip(cos_d, -1, 1)
    # Earth's radius in meters
    R = 6371000
    distance = R * np.arccos(cos_d)

    # Set distance to zero for identical points
    identical = np.all(np.abs(latlon1[:, np.newaxis] - latlon2) < 1e-10, axis=2)
    distance[identical] = 0

    return distance


class GeospatialLocator(object):
    """
    Class to interact with a postcode database file and a population grid file.
    """

    def __init__(
        self,
        postcode_file=os.sep.join(
            (os.path.dirname(__file__), "..", "resources", "full_postcodes.csv")
        ),
        census_file=os.sep.join(
            (
                os.path.dirname(__file__),
                "..",
                "resources",
                "UK_residential_population_2011_latlon.asc",
            )
        ),
        norm=great_circle_distance,
    ):
        """
        Parameters
        ----------

        postcode_file : str, optional
            Filename of a .csv file containing geographic
            location data for postcodes.

        census_file :  str, optional
            Filename of a .asc file containing census data on a
            latitude-longitude grid.

        norm : function
            Python function defining the distance between points in
            latitude-longitude space.

        """

        self.postcode_file = postcode_file
        self.census_file = census_file
        self.norm = norm

        with open(self.census_file, "r") as file:
            # Headers
            ncols, nrows = int(file.readline().split()[1]), int(
                file.readline().split()[1]
            )
            nodata_value = float(file.readline().split()[1])
            # Skipping the array labels
            for _ in range(3):
                file.readline()

            # Read the data
            data = np.loadtxt(file, dtype=float)
            # Reshape into three stacked arrays: latitude, longitude, population
            data = data.reshape((3, nrows, ncols))

            # Updating the population to 0 for missing values
            data[2][data[2] == nodata_value] = 0

            latitude, longitude, population = [arr.flatten() for arr in data]
        self.census = pd.DataFrame(
            {"Latitude": latitude, "Longitude": longitude, "Population": population}
        )

    def get_postcodes_by_radius(self, X, radii):
        """
        Return postcodes within specific distances of
        input location.

        Parameters
        ----------
        X : arraylike
            Latitude-longitude pair of centre location
        radii : arraylike
            array of radial distances from X

        Returns
        -------
        list of lists
            Contains the lists of postcodes closer than the elements
            of radii to the location X.


        Examples
        --------

        >>> locator = GeospatialLocator()
        >>> locator.get_postcodes_by_radius((51.4981, -0.1773), [1.5e3])
        >>> locator.get_postcodes_by_radius((51.4981, -0.1773),
                                            [1.5e3, 4.0e3])
        """
        # Load postcode data and filter valid coordinates
        postcodes = pd.read_csv(self.postcode_file).query(
            "-90 <= Latitude <= 90 and -180 <= Longitude <= 180"
        )
        # Return an empty list for empty postcodes
        if postcodes.empty:
            return [[] for _ in radii]

        # Calculating distances to all postcodes
        distances = self.norm(postcodes[["Latitude", "Longitude"]].values, [X])

        result = []
        for radius in radii:
            if radius <= 0:
                # Return an empty list for non-positive radius values
                result.append([])
                continue

            # Filter postcodes within the radius
            within_radius = postcodes[distances.squeeze() <= radius]
            result.append(list(within_radius["Postcode"]))

        return result

    def find_nearest_coords(self, X, num_coords):
        """
        Find the nearest census grid points to a given location.

        Parameters
        ----------
        X : arraylike
            Latitude-longitude pair [lat, lon] of the reference point.
        num_coords : int
            Number of nearest coordinates to find.

        Returns
        -------
        coords : ndarray
            Array of shape (num_coords, 2) containing nearest coordinates.
        dist : list
            Great circle distances (in metres) to each coordinate.

        Notes
        -----
        Uses spatial filtering with expanding search radius for efficient
        nearest neighbor search on 1km×1km grid data.

        The method progressively expands the search area until enough points
        are found, then performs brute-force distance calculation within that area.

        Time complexity: O(k) on average for uniformly distributed grid data,
        where k is proportional to the number of points in the search area.
        """
        lat_query, lon_query = X

        # Grid is 1km×1km, convert km to approximate degrees
        km_to_lat = 0.009  # degrees per km
        km_to_lon = 0.014 / np.cos(np.radians(lat_query))  # adjust for latitude

        # Start with a small search radius and expand if needed
        # Initial radius: enough to cover sqrt(num_coords) × sqrt(num_coords) grid
        search_radius_km = max(3, np.ceil(np.sqrt(num_coords)) + 2)
        max_iterations = 10

        for iteration in range(max_iterations):
            lat_margin = search_radius_km * km_to_lat
            lon_margin = search_radius_km * km_to_lon

            lat_min = lat_query - lat_margin
            lat_max = lat_query + lat_margin
            lon_min = lon_query - lon_margin
            lon_max = lon_query + lon_margin

            # Filter census data within bounding box
            candidates = self.census[
                (self.census["Latitude"] >= lat_min)
                & (self.census["Latitude"] <= lat_max)
                & (self.census["Longitude"] >= lon_min)
                & (self.census["Longitude"] <= lon_max)
            ]

            # If we have enough candidates, compute distances and find nearest
            if len(candidates) >= num_coords:
                # Brute force distance calculation within the filtered area
                candidate_coords = candidates[["Latitude", "Longitude"]].values
                distances = self.norm([X], candidate_coords)[0]

                nearest_indices = np.argpartition(
                    distances, min(num_coords - 1, len(distances) - 1)
                )[:num_coords]
                # Sort these indices by distance
                nearest_indices = nearest_indices[
                    np.argsort(distances[nearest_indices])
                ]

                coords = candidate_coords[nearest_indices]
                dist = distances[nearest_indices].tolist()

                return coords, dist

            # Not enough points, expand search radius
            search_radius_km *= 2

        candidate_coords = candidates[["Latitude", "Longitude"]].values
        distances = self.norm([X], candidate_coords)[0]

        # Sort all candidates by distance
        sorted_indices = np.argsort(distances)[:num_coords]
        coords = candidate_coords[sorted_indices]
        dist = distances[sorted_indices].tolist()

        return coords, dist

    def cal_impacted_population(self, radius, grid_centers, grid_distances):
        """
        Calculate the population impacted within a specified radius around a set of grid centers.

        Parameters
        ----------
        radius : float
            The radius within which to calculate the impacted population.
        grid_centers : arraylike
            An array of latitude-longitude pairs representing grid centers.
        grid_distances : arraylike
            An array of distances from a reference point to each grid center.

        Returns
        -------
        float
            The total population impacted within the specified radius around the given grid centers.

        Notes
        -----
        This method approximates the impacted population by considering the proportion of each grid cell's area
        that falls within the specified radius. The population data is sourced from the census data loaded earlier.
        """
        # Initialize the impacted population
        impacted_pop = 0
        # Calculate the half-diagonal of the square grid
        half_diagonal = 1000 / np.sqrt(2)

        for grid_center, distance in zip(grid_centers, grid_distances):
            intersection_percentage = 0
            population = self.census[
                (self.census["Latitude"] == grid_center[0])
                & (self.census["Longitude"] == grid_center[1])
            ]["Population"].values[0]

            # Check if the entire circle is within a grid cell
            if distance + radius < half_diagonal:
                intersection_percentage = (np.pi * radius**2) / (1000 * 1000)
                impacted_pop += intersection_percentage * population
                return impacted_pop

            # Check if the grid intersects with the circle
            if distance <= radius + half_diagonal:
                # Calculate intersection percentage
                if distance <= radius - half_diagonal:
                    intersection_percentage = 1
                else:
                    intersection_percentage = max(
                        0,
                        1 - (distance - (radius - half_diagonal)) / (2 * half_diagonal),
                    )

                impacted_pop += intersection_percentage * population

        return impacted_pop

    def get_population_by_radius(self, X, radii, method="naive"):
        """
        Return the population within specific distances of input location.

        Parameters
        ----------
        X : arraylike
            Latitude-longitude pair of centre location
        radii : arraylike
            array of radial distances from X

        Returns
        -------
        list
            Contains the population closer than the elements of radii to
            the location X. Output should be the same shape as the radii array.

        Examples
        --------
        >>> loc = GeospatialLocator()
        >>> loc.get_population_by_radius((51.4981, -0.1773), [1e2, 5e2, 1e3])

        """
        # Calculate distances from X to each point in the census data
        self.census["distance"] = self.norm(
            self.census[["Latitude", "Longitude"]].values, [X]
        )[:, 0]

        populations_by_radius = []
        for radius in radii:
            if radius <= 0:
                populations_by_radius.append(0)
                continue

            # Sum population for points within the radius
            if radius <= 500:
                # Find the 4 nearest coordinates
                coords, dist = self.find_nearest_coords(X, 4)
                total_population = self.cal_impacted_population(radius, coords, dist)

            elif radius > 500 and radius <= 1000:
                # Find the 10 nearest coordinates
                coords, dist = self.find_nearest_coords(X, 10)
                total_population = self.cal_impacted_population(radius, coords, dist)

            else:
                total_population = self.census[self.census["distance"] <= radius][
                    "Population"
                ].sum()

            populations_by_radius.append(int(total_population))

        return populations_by_radius
