from slice import generate_slices
from fileio import load_globals
from plot import plot_slices

data = load_globals("input_template.xlsx")
profile_lines = data["profile_lines"]
materials = data["materials"]
piezo_line = data["piezo_line"]
circle = data["circles"][0]  # or whichever one you want

df, failure_surface = generate_slices(profile_lines, materials, circle, piezo_line=piezo_line)
print(df)
print(failure_surface)

# plot_slices(profile_lines, circle, df, piezo_line=piezo_line, failure_surface=failure_surface)
