from xslope.fileio import load_slope_data, print_dictionary

# A reinforced model, so the dump covers materials, profile lines, polygons,
# piezometric data and reinforcement lines rather than just the geometry.
filepath = "docs/inputs/slope/xslope_reinf.xlsx"
slope_data = load_slope_data(filepath)

print_dictionary(slope_data)
