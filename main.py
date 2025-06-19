from global_config import non_circ
from slice import generate_slices
from fileio import load_globals
from plot import plot_solution, plot_inputs
from solve import oms, bishop, janbu, spencer, corps_engineers, lowe_karafiath


def solve_selected(func, df):
    success, result = func(df)
    if not success:
        print(f'Error: {result}')
        return result

    if func == oms:
        print(f'OMS: FS={result["FS"]:.3f}')
    elif func == bishop:
        print(f'Bishop: FS={result["FS"]:.3f}')
    elif func == spencer:
        print(f'Spencer: FS={result["FS"]:.3f}, theta={result["theta"]:.2f}')
    elif func == janbu:
        print(f'Janbu Corrected FS={result["FS"]:.3f}, fo={result["fo"]:.2f}')
    elif func == corps_engineers:
        print(f'Corps Engineers: FS={result["FS"]:.3f}, theta={result["theta"]:.2f}')
    elif func == lowe_karafiath:
        print(f'Lowe & Karafiath: FS={result["FS"]:.3f}')
    return result

def solve_all(df):
    solve_selected(oms, df)
    solve_selected(bishop, df)
    solve_selected(janbu, df)
    solve_selected(corps_engineers, df)
    solve_selected(lowe_karafiath, df)
    solve_selected(spencer, df)

data = load_globals("docs/input_template_dam.xlsx")

# plot_inputs(data)

circle = data['circles'][0] if data['circular'] else None
# non_circ = data['non_circ'] if data['non_circ'] else None

# print(f"circle: {circle}")

success, result = generate_slices(data, circle=circle, non_circ=None, num_slices=20)
# success, result = generate_slices(data, circle=None, non_circ=non_circ, num_slices=20)
if success:
    df, failure_surface = result
else:
    print(result)

# export df to excel
df.to_csv("slices.csv", index=False)

# options = [oms, bishop, janbu, corps_engineers, lowe_karafiath, spencer]
# results = solve_selected(oms, df)

solve_all(df)

# plot_solution(data, df, failure_surface, results)

