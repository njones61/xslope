from fileio import load_globals

from mesh import build_polygons, print_polygon_summary, plot_polygons
from plot import plot_inputs



data = load_globals("inputs/slopes/input_template_lface2.xlsx")

plot_inputs(data)

polygons = build_polygons(data['profile_lines'], max_depth=data['max_depth'])

print_polygon_summary(polygons)

plot_polygons(polygons)


