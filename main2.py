from slice import generate_slices
from fileio import load_globals
from plot import plot_slices
from utils import ground_surface

data = load_globals("input_template.xlsx")
profile_lines = data["profile_lines"]
materials = data["materials"]
piezo_line = data["piezo_line"]
circle = data["circles"][0]  # or whichever one you want

top_surface = ground_surface(profile_lines)

print(top_surface)

df, arc = generate_slices(profile_lines, materials, circle, piezo_line, num_slices=20, surface_polyline=top_surface)


print(df[df.columns[:15]])  # Adjust the slicing as needed

plot_slices(profile_lines, circle, df, piezo_line=piezo_line, failure_surface=arc)
