from fileio import load_slope_data
from plot import plot_inputs

# Test with a file that has reinforcement data
filepath = "inputs/slope/input_template_reinf5.xlsx"
slope_data = load_slope_data(filepath)

# Test the updated plotting function
print("Testing updated plot_inputs function with smaller markers and legend...")
plot_inputs(slope_data, title="Test Plot with Updated Reinforcement Lines", mat_table=True) 