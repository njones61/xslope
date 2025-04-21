##################################################################
#### This file contains global variables used in the project. ####
##################################################################

"""Global variables used throughout the project, documented for mkdocstrings."""

#: [lb/ft^3] Unit weight of water
gamma_water = 62.4

#: [ft] Depth of the crack
tcrack_depth = 0.0

#: [ft] Water level in the crack
tcrack_water = 0.0

#: List of profile lines, each a list of (x, y) tuples.
profile_lines = []
# Example: profile_lines = [[(x1, y1), (x2, y2)], [(x3, y3), (x4, y4)]]

#: List of material property dictionaries.
materials = []
# Example: materials = [{'gamma': 120, 'c': 30, 'phi': 20, 'piezo': 0.5, 'sigma_gamma': 0.1, 'sigma_c': 0.1, 'sigma_phi': 0.1}]

#: List of (x, y) coordinates defining the piezometric line.
piezo_line = []
# Example: piezo_line = [(x1, y1), (x2, y2)]

#: [ft] Maximum depth for circular failure surfaces.
max_depth = 0.0

#: List of circular failure surface dictionaries.
circles = []
# Example: circles = [{'Xo': 120, 'Yo': 80, 'Option': "Depth", 'Depth': -10, 'Xi': 5, 'Yi': 5}]

#: List of points in a non-circular surface with movement options.
non_circ = []
# Example: non_circ = [{'X': 120, 'Y': 80, 'Movement': "Free"}, {'X': 130, 'Y': 90, 'Movement': "Horiz"}]

#: List of distributed load dictionaries.
dloads = []
# Example: dloads = [{'X': 120, 'Y': 80, 'Normal': 100}, {'X': 130, 'Y': 90, 'Normal': 150}]

#: List of reinforcement line dictionaries.
reinforce_lines = []
# Example: reinforce_lines = [{'X': 120, 'Y': 80, 'FL': 10, 'FT': 20}, {'X': 130, 'Y': 90, 'FL': 15, 'FT': 25}]
