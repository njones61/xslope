from fileio import load_data_from_pickle, save_data_to_pickle, load_globals
import os

filepath = "docs/input_template_lface2"  # without extensions

print("Loading data from Excel file...")
data = load_globals(filepath + ".xlsx")
print(f"Data loaded successfully. Keys: {list(data.keys())}")

print("Saving data to pickle file...")
save_data_to_pickle(data, filepath + ".pkl")
print("Pickle file saved successfully!")

print("Verifying file was created...")
if os.path.exists(filepath + ".pkl"):
    print(f"File exists! Size: {os.path.getsize(filepath + '.pkl')} bytes")
else:
    print("File was not created!")







