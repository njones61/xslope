
from fileio import load_slope_data

filepath = "inputs/slope/input_template_MASTER5.xlsx"  # Replace with full path if needed
slope_data = load_slope_data(filepath)

for key, value in slope_data.items():
    print(f"\n=== {key} ===")
    if isinstance(value, list):
        for item in value:
            print(item)
    else:
        print(value)
