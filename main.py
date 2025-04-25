from global_config import non_circ
from slice import generate_slices
from fileio import load_globals
from plot import plot_slope
from utils import build_ground_surface
import math
from solve import oms, bishop, spencer, janbu_corrected


def solve_selected(func, df, circular=True):
    results = func(df, circular=circular)
    if func == oms:
        print(f'OMS: FS={results["FS"]:.3f}')
    elif func == bishop:
        print(f'Bishop: FS={results["FS"]:.3f}')
    elif func == spencer:
        print(f'Spencer: FS={results["FS"]:.3f}, theta={results["theta"]:.2f}')
    elif func == janbu_corrected:
        print(f'Janbu Corrected FS={results["FS"]:.3f}, fo={results["fo"]:.2f}')
    if not results['success']:
        print(f'Error: {results["message"]}')
    return results

def solve_all(df, circular=True):
    solve_selected(oms, df, circular=circular)
    solve_selected(bishop, df, circular=circular)
    solve_selected(oms, df, circular=circular)
    solve_selected(spencer, df, circular=circular)
    solve_selected(janbu_corrected, df, circular=circular)


data = load_globals("docs/input_template.xlsx")

# plot_slope(data)

ground_surface = build_ground_surface(data['profile_lines'])
circle = data['circles'][0] if data['circular'] else None
non_circ = data['non_circ'] if not data['circular'] else None

success, result = generate_slices(data, ground_surface, circle, non_circ, num_slices=20)
if success:
    df, failure_surface = result
else:
    print(result)

# export df to excel
#df.to_excel("slices.xlsx", index=False)


# options = [oms, bishop, spencer, janbu_corrected]
results = solve_selected(oms, df, circular=True)
plot_slope(data, df=df, failure_surface=failure_surface, fs=results['FS'])

# solve_all(df, circular=True)