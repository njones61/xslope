from slice import generate_slices
from fileio import load_globals

data = load_globals("input_template.xlsx")
profile_lines = data["profile_lines"]
materials = data["materials"]
circle = data["circles"][0]  # or whichever one you want

df = generate_slices(profile_lines, materials, circle)
print(df)