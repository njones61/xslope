from global_config import non_circ
from slice import generate_slices
from fileio import load_slope_data, load_data_from_pickle
from plot import plot_solution, plot_inputs
from solve import oms, bishop, janbu, corps_engineers, lowe_karafiath, spencer, spencer_OLD


def solve_selected(func, slice_df, rapid=False):
    if rapid:
        success, result = rapid_drawdown(slice_df, func)
    else:
        success, result = func(slice_df)
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

def solve_all(slice_df):
    solve_selected(oms, slice_df)
    solve_selected(bishop, slice_df)
    solve_selected(janbu, slice_df)
    solve_selected(corps_engineers, slice_df)
    solve_selected(lowe_karafiath, slice_df)
    solve_selected(spencer, slice_df)

slope_data = load_slope_data("inputs/slope/input_template_reliability4.xlsx")

# plot_inputs(slope_data)

circle = slope_data['circles'][0] if slope_data['circular'] else None
non_circ = slope_data['non_circ'] if slope_data['non_circ'] else None

success, result = generate_slices(slope_data, circle=circle, non_circ=None, num_slices=20)

if success:
    slice_df, failure_surface = result
else:
    print(result)
    exit()

# options = [oms, bishop, janbu, corps_engineers, lowe_karafiath, spencer]
results = solve_selected(spencer, slice_df, rapid=False)

# solve_all(slice_df)

# slice_df.to_excel("slices.xlsx", index=False)

plot_solution(slope_data, slice_df, failure_surface, results, slice_numbers=True)

