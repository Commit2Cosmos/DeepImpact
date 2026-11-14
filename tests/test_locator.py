import numpy as np

from pytest import fixture, mark


# Use pytest fixtures to generate objects we know we'll reuse.
# This makes sure tests run quickly


@fixture(scope="module")
def deepimpact():
    import deepimpact

    return deepimpact


@fixture(scope="module")
def loc(deepimpact):
    return deepimpact.GeospatialLocator()


@mark.xfail
def test_great_circle_distance(deepimpact):

    pnts1 = np.array([[54.0, 0.0], [55.0, 1.0], [54.2, -3.0]])
    pnts2 = np.array([[55.0, 1.0], [56.0, -2.1], [54.001, -0.003]])

    data = np.array(
        [
            [1.28580537e05, 2.59579735e05, 2.25409117e02],
            [0.00000000e00, 2.24656571e05, 1.28581437e05],
            [2.72529953e05, 2.08175028e05, 1.96640630e05],
        ]
    )

    dist = deepimpact.great_circle_distance(pnts1, pnts2)

    assert np.allclose(data, dist, rtol=1.0e-4)


def test_locator_postcodes(loc):

    latlon = (58.2074, 8.1170)

    result = loc.get_postcodes_by_radius(latlon, [1.2e3, 1.1e3, 1.3e5])

    assert type(result) is list
    if len(result) > 0:
        for element in result:
            assert type(element) is list


def test_population_by_radius(loc):

    latlon = (52.2074, 0.1170)

    result = loc.get_population_by_radius(latlon, [4e2, 5e7, 1e10], method="naive")

    assert type(result) is list
    if len(result) > 0:
        for element in result:
            assert type(element) is int


def test_find_nearest_coords(loc):
    """Test finding nearest census grid coordinates."""

    latlon = (52.2074, 0.1170)
    num_coords = 5

    coords, distances = loc.find_nearest_coords(latlon, num_coords)

    # Check return types
    assert isinstance(coords, np.ndarray)
    assert isinstance(distances, list)

    # Check dimensions
    assert coords.shape == (num_coords, 2)
    assert len(distances) == num_coords


def test_cal_impacted_population(loc):
    """Test calculating impacted population for given grid centers."""

    radius = 400

    # Use find_nearest_coords to get actual grid centers
    latlon = (52.2074, 0.1170)
    grid_centers, grid_distances = loc.find_nearest_coords(latlon, 4)

    result = loc.cal_impacted_population(radius, grid_centers, grid_distances)

    # Test with large radius should return more population
    result_large = loc.cal_impacted_population(800, grid_centers, grid_distances)
    assert result_large >= result


def test_population_by_radius_ranges(loc):
    """Test that population increases with radius."""

    latlon = (52.2074, 0.1170)

    # Test small radius (< 500m)
    result_small = loc.get_population_by_radius(latlon, [200, 400])
    assert all(isinstance(p, int) for p in result_small)
    assert result_small[1] >= result_small[0]  # Larger radius >= smaller radius

    # Test medium radius (500-1000m)
    result_medium = loc.get_population_by_radius(latlon, [600, 800])
    assert all(isinstance(p, int) for p in result_medium)
    assert result_medium[1] >= result_medium[0]

    # Test large radius (> 1000m)
    result_large = loc.get_population_by_radius(latlon, [2000, 5000])
    assert all(isinstance(p, int) for p in result_large)
    assert result_large[1] >= result_large[0]
