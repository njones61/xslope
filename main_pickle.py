from fileio import load_data_from_pickle, save_data_to_pickle, load_slope_data
import os

filepath = "docs/input_template_lface2"  # without extensions

print("Loading data from Excel file...")
slope_data = load_slope_data(filepath + ".xlsx")
print(f"Data loaded successfully. Keys: {list(slope_data.keys())}")

print("Saving data to pickle file...")
save_data_to_pickle(slope_data, filepath + ".pkl")
print("Pickle file saved successfully!")

print("Verifying file was created...")
if os.path.exists(filepath + ".pkl"):
    print(f"File exists! Size: {os.path.getsize(filepath + '.pkl')} bytes")
else:
    print("File was not created!")







