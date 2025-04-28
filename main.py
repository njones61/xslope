from global_config import non_circ
from slice import build_ground_surface, generate_slices
from fileio import load_globals
from plot import plot_slope
from solve import oms, bishop, spencer, janbu_corrected


def solve_selected(func, df, circular=True):
    success, result = func(df, circular=circular)
    if not success:
        print(f'Error: {result}')
        return result

    if func == oms:
        print(f'OMS: FS={result["FS"]:.3f}')
    elif func == bishop:
        print(f'Bishop: FS={result["FS"]:.3f}')
    elif func == spencer:
        print(f'Spencer: FS={result["FS"]:.3f}, theta={result["theta"]:.2f}')
    elif func == janbu_corrected:
        print(f'Janbu Corrected FS={result["FS"]:.3f}, fo={result["fo"]:.2f}')
    return result

def solve_all(df, circular=True):
    solve_selected(oms, df, circular=circular)
    solve_selected(bishop, df, circular=circular)
    solve_selected(oms, df, circular=circular)
    solve_selected(spencer, df, circular=circular)
    solve_selected(janbu_corrected, df, circular=circular)


data = load_globals("docs/input_template.xlsx")

# plot_slope(data)

circle = data['circles'][0] if data['circular'] else None
non_circ = data['non_circ'] if not data['circular'] else None

success, result = generate_slices(data, circle, non_circ, num_slices=20)
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