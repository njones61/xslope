from slice import generate_slices
from fileio import load_globals
from plot import plot_slices
from utils import ground_surface
from solve import oms, bishops, spencer, janbu_simple, janbu_corrected, morgenstern_price

data = load_globals("input_template.xlsx")
profile_lines = data["profile_lines"]
materials = data["materials"]
piezo_line = data["piezo_line"]
gamma_w = data["gamma_water"]
circle = data["circles"][0]  # or whichever one you want
dload = data["dloads"]

top_surface = ground_surface(profile_lines)
df, arc = generate_slices(profile_lines, materials, circle,  top_surface, num_slices=20, gamma_w=gamma_w, piezo_line=piezo_line, dloads=dload)

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

# FS, beta_deg = spencer(df)
# print(f"Spencer's Method: FS = {FS:.4f}, β = {beta_deg:.2f}°")

# FS = janbu_simple(df)
# print(f"Janbu FS = {FS:.4f}")

# FS, lam, converged = janbu_corrected(df)
# if converged:
#     print(f"Janbu Corrected FS = {FS:.4f}, λ = {lam:.4f}")
# else:
#     print(f"Did not converge. Last FS = {FS:.4f}, λ = {lam:.4f}")




# Default constant interslice force function (Janbu)
FS, lam, converged = morgenstern_price(df)

# With a half-sine function:
# import math
# FS, lam, converged = morgenstern_price(df, function=lambda x: math.sin(math.pi * x))

if converged:
    print(f"Morgenstern-Price FS = {FS:.4f}, λ = {lam:.4f}")
else:
    print(f"Did not converge. Last FS = {FS:.4f}, λ = {lam:.4f}")

plot_slices(profile_lines, circle, df, piezo_line=piezo_line, failure_surface=arc, FS=FS, dloads=dload)
