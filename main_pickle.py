import os
import tempfile

from xslope.fileio import load_slope_data, save_data_to_pickle, load_data_from_pickle

xlsx_path = "docs/inputs/slope/xslope_lface.xlsx"

print("Loading data from Excel file...")
slope_data = load_slope_data(xlsx_path)
print(f"Data loaded successfully. Keys: {list(slope_data.keys())}")

# A pickle round trip is faster than re-parsing the workbook, which matters when
# the same model is solved repeatedly. Written to a temp file so the demo leaves
# nothing behind; in real use it sits beside the .xlsx.
pkl_path = os.path.join(tempfile.gettempdir(), "xslope_lface.pkl")

print("Saving data to pickle file...")
save_data_to_pickle(slope_data, pkl_path)
print(f"Pickle file saved: {os.path.getsize(pkl_path)} bytes")

print("Reading it back...")
reloaded = load_data_from_pickle(pkl_path)
print(f"Reloaded keys match: {list(reloaded.keys()) == list(slope_data.keys())}")

os.remove(pkl_path)
