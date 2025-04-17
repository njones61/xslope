##################################################################
#### This file contains global variables used in the project. ####
##################################################################

gamma_water = 62.4    # [lb/ft^3] Unit weight of water
tcrack_depth = 0.0    # [ft] Depth of the crack
tcrack_water = 0.0    # [ft] Water level in the crack

profile_lines = []    # List to store profile lines. Each line is a list of XY coordinates (tuples).
# Each profile line is a list of tuples, where each tuple represents a point (x, y).
# Example: profile_lines = [[(x1, y1), (x2, y2)], [(x3, y3), (x4, y4)]]

materials = []    # List to store material properties. Each material is a dictionary with keys:
# gamma, option, c, phi, cp, r_elev piezo, sigma_gamma, sigma_c, sigma_phi, sigma_cp.
# Example: materials = [{'gamma': 120, 'c': 30, 'phi': 20, 'piezo': 0.5, 'sigma_gamma': 0.1, 'sigma_c': 0.1, 'sigma_phi': 0.1},

piezo_line = []    # List to store piezometric line defined by a list of XY coordinates (tuples).
# Example: piezo_line = [(x1, y1), (x2, y2)]

max_depth = 0.0    # [ft] Maximum depth that circles can reach

circles = []    # List to store circles. Each circle is a dictionary with keys:
# Xo, Yo, Option, Depth, Xi, Yi
# where Option is either "Depth" or "Intercept".
# Example: circles = [{'Xo': 120, 'Yo': 80, 'Option': "Depth", 'Depth': -10, 'Xi': 5, 'Yi': 5}

non_circ = []    # List to store coordinates of a non-circular shape and a code for each point.
# Each non-circular shape is a dictionary with keys:
# X, Y, Movement
# The options for movement are: "Free" and "Horiz".
# Example: non_circ = [{'X': 120, 'Y': 80, 'Movement': "Free"}, {'X': 130, 'Y': 90, 'Movement': "Horiz"}]

dloads = []    # List to store distributed loads. Each load is a dictionary with keys:
# X, Y, Normal
# Example: dloads = [{'X': 120, 'Y': 80, 'Normal': 100}, {'X': 130, 'Y': 90, 'Normal': 150}]

reinforce_lines = [] # List to store reinforcement lines. Each line is a dictionary with keys:
# X, Y, FL, FT
# Example: reinforce_lines = [{'X': 120, 'Y': 80, 'FL': 10, 'FT': 20}, {'X': 130, 'Y': 90, 'FL': 15, 'FT': 25}]
# where FL is the longitudinal forces and FT is the transverse force.

