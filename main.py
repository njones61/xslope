from global_config import non_circ
from slice import generate_slices
from fileio import load_globals
from plot import plot_slices
from utils import build_ground_surface
from solve import oms, bishop, spencer, spencer_moment, janbu_simple, janbu_corrected, morgenstern_price


def solve_selected(method):
    # Options: 'oms', 'bishop', 'spencer',  spencer_moment, 'janbu_simple', 'janbu_corrected', 'morgenstern_price (janbu)', 'morgenstern_price (half-sine)'

    if method == 'oms':
        FS, N = oms(df)
        print(f"OMS: FS = {FS:.3f}")
    elif method == 'bishop':
        FS, N, converge = bishop(df)
        if not converge:
            print("Bishop's method did not converge.")
        print(f"Bishop: FS = {FS:.3f}")
    elif method == 'spencer':
        FS, beta_deg = spencer(df)
        print(f"Spencer's Method: FS = {FS:.4f}, β = {beta_deg:.2f}°")
    elif method == 'spencer_moment':
        FS, beta_deg = spencer_moment(df)
        print(f"Spencer's Method (moment): FS = {FS:.4f}, β = {beta_deg:.2f}°")
    elif method == 'janbu_simple':
        FS = janbu_simple(df)
        print(f"Janbu's Method: FS = {FS:.4f}")
    elif method == 'janbu_corrected':
        FS, lam, converged = janbu_corrected(df)
        if converged:
            print(f"Janbu's Corrected Method: FS = {FS:.4f}, λ = {lam:.4f}")
        else:
            print(f"Did not converge. Last FS = {FS:.4f}, λ = {lam:.4f}")
    elif method == 'morgenstern_price (janbu)':
        # Default constant interslice force function (Janbu)
        FS, lam, converged = morgenstern_price(df)
        if converged:
            print(f"Morgenstern-Price Method: FS = {FS:.4f}, λ = {lam:.4f}")
        else:
            print(f"Did not converge. Last FS = {FS:.4f}, λ = {lam:.4f}")
    elif method == 'morgenstern_price (half-sine)':
        # With a half-sine function:
        import math
        FS, lam, converged = morgenstern_price(df, function=lambda x: math.sin(math.pi * x))
        if converged:
            print(f"Morgenstern-Price Method (half-sine): FS = {FS:.4f}, λ = {lam:.4f}")
        else:
            print(f"Did not converge. Last FS = {FS:.4f}, λ = {lam:.4f}")
    return FS

data = load_globals("docs/input_template.xlsx")
profile_lines = data["profile_lines"]
materials = data["materials"]
piezo_line = data["piezo_line"]
gamma_w = data["gamma_water"]
circle = data["circles"][0]  # or whichever one you want
non_circ = data["non_circ"]
dloads = data["dloads"]
max_depth = data["max_depth"]
reinforce_lines = data["reinforce_lines"]

ground_surface = build_ground_surface(profile_lines)

df, failure_surface = generate_slices(
    profile_lines=profile_lines,
    materials=materials,
    ground_surface=ground_surface,
    circle=circle,
    #non_circ=non_circ,
    num_slices=20,
    gamma_w=62.4,
    piezo_line=piezo_line,
    dloads=dloads,
    #reinforce_lines=reinforce_lines
)


# export df to excel
df.to_excel("slices.xlsx", index=False)

#print(df[df.columns[10:]])

#print(df[df.columns[0,]])  # Adjust the slicing as needed

# options = ['oms', 'bishop', 'spencer',  'spencer_moment', 'janbu_simple', 'janbu_corrected', 'morgenstern_price (janbu)', 'morgenstern_price (half-sine)']
method = 'spencer'  # Change this to the desired method
FS = solve_selected(method)


plot_slices(profile_lines, df, piezo_line=piezo_line, failure_surface=failure_surface, fs=FS, dloads=dloads, max_depth=max_depth)
