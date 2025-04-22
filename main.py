from global_config import non_circ
from slice import generate_slices
from fileio import load_globals
from plot import plot_slices
from utils import build_ground_surface
import math
from solve import oms, bishop, spencer, janbu_corrected


def solve_selected(method, circular=False):

    if method == 'oms':
        FS, N = oms(df)
        print(f"OMS: FS = {FS:.3f}")
    elif method == 'bishop':
        FS, N, converge = bishop(df)
        if not converge:
            print("Bishop's method did not converge.")
        print(f"Bishop: FS = {FS:.3f}")
    elif method == 'spencer':
        FS, beta_deg, converged = spencer(df, circular=circular)
        if not converged:
            print("Spencer's method did not converge.")
        print(f"Spencer's Method (SGW): FS = {FS:.4f}, β = {beta_deg:.2f}°")
    elif method == 'janbu_corrected':
        FS, fo, converged = janbu_corrected(df)
        if not converged:
            print("Janbu's corrected method did not converge.")
        else:
            print(f"Janbu's Corrected Method: FS = {FS:.4f}, fo = {fo:.4f}")
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
    #circle=circle,
    non_circ=non_circ,
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

# options = ['oms', 'bishop', 'spencer', 'janbu_corrected']
method = 'janbu_corrected'  # Change this to the desired method
FS = solve_selected(method, circular=False)


plot_slices(profile_lines, df, piezo_line=piezo_line, failure_surface=failure_surface, fs=FS, dloads=dloads, max_depth=max_depth)
