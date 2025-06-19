from fileio import load_globals
from shapely.geometry import LineString
from slice import get_y_from_intersection

data = load_globals('docs/input_template_dam.xlsx')
print('Profile line intersections at different x-coordinates:')

print(data['ground_surface'])