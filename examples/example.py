import deepimpact

#######################
#   Airburst Solver   #
#######################

# Initialise the Planet class
earth = deepimpact.Planet()

# Solve the atmospheric entry problem for a given set of input parameters
result = earth.solve_atmospheric_entry(radius=35, angle=45,
                                       strength=1e7, density=3000,
                                       velocity=19e3)

# Calculate the kinetic energy lost per unit altitude and add it
# as a column to the result dataframe
result = earth.calculate_energy(result)

# Determine the outcomes of the impact event
outcome = earth.analyse_outcome(result)
