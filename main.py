from slice import generate_slices
from fileio import load_globals
from plot import plot_slices
from utils import ground_surface
from solve import oms, bishops, spencer

data = load_globals("input_template.xlsx")
profile_lines = data["profile_lines"]
materials = data["materials"]
piezo_line = data["piezo_line"]
gamma_w = data["gamma_water"]
circle = data["circles"][0]  # or whichever one you want

top_surface = ground_surface(profile_lines)
df, arc = generate_slices(profile_lines, materials, circle,  top_surface, num_slices=20, gamma_w=gamma_w, piezo_line=piezo_line)

# export df to excel
df.to_excel("slices.xlsx", index=False)

#print(df[df.columns[10:]])

#print(df[df.columns[0,]])  # Adjust the slicing as needed

# FS, N = oms(df)
# print(f"Factor of Safety = {FS:.3f}")


# FS, N, converge = bishops(df)
# if not converge:
#     print("Bishop's method did not converge.")
# print(f"Factor of Safety = {FS:.3f}")

FS, beta_deg = spencer(df)
print(f"Spencer's Method: FS = {FS:.4f}, β = {beta_deg:.2f}°")

plot_slices(profile_lines, circle, df, piezo_line=piezo_line, failure_surface=arc, FS=FS)
