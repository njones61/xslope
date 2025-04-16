
import pandas as pd

def load_globals(filepath):
    xls = pd.ExcelFile(filepath)
    globals_data = {}

    # === PROFILE LINES ===
    profile_df = xls.parse('profile', header=None)
    profile_lines = []

    profile_data_blocks = [
        {"header_row": 2, "data_start": 3, "data_end": 18},   # Set 1
        {"header_row": 20, "data_start": 21, "data_end": 36}  # Set 2
    ]
    profile_block_width = 3

    for block in profile_data_blocks:
        for col in range(0, profile_df.shape[1], profile_block_width):
            x_col, y_col = col, col + 1
            try:
                x_header = str(profile_df.iloc[block["header_row"], x_col]).strip().lower()
                y_header = str(profile_df.iloc[block["header_row"], y_col]).strip().lower()
            except:
                continue
            if x_header != 'x' or y_header != 'y':
                continue
            data = profile_df.iloc[block["data_start"]:block["data_end"], [x_col, y_col]].dropna()
            coords = data.apply(lambda r: (float(r.iloc[0]), float(r.iloc[1])), axis=1).tolist()
            if coords:
                profile_lines.append(coords)

    globals_data["profile_lines"] = profile_lines

    # === MATERIALS ===
    mat_df = xls.parse('mat')
    materials = []
    for _, row in mat_df.iloc[2:].iterrows():
        try:
            gamma, c, phi = float(row[1]), float(row[2]), float(row[3])
        except:
            continue
        materials.append({
            "gamma": gamma,
            "c": c,
            "phi": phi,
            "piezo": float(row[4]),
            "sigma_gamma": float(row[6]) if pd.notna(row[6]) else 0,
            "sigma_c": float(row[7]) if pd.notna(row[7]) else 0,
            "sigma_phi": float(row[8]) if pd.notna(row[8]) else 0,
        })
    globals_data["materials"] = materials

    # === PIEZOMETRIC LINE ===
    piezo_df = xls.parse('piezo')
    piezo_line = list(piezo_df.iloc[2:].dropna().apply(lambda row: (float(row[0]), float(row[1])), axis=1))
    globals_data["piezo_line"] = piezo_line

    # === CIRCLES ===
    circles_df = xls.parse('circles')
    globals_data["max_depth"] = float(circles_df.iloc[0, 2])
    circles_rows = circles_df.iloc[3:]
    circles_rows = circles_rows.dropna(subset=circles_df.columns[[1, 2]])
    circles = []
    for _, row in circles_rows.iterrows():
        circles.append({
            "Xo": float(row[1]),
            "Yo": float(row[2]),
            "Option": row[3],
            "Depth": float(row[4]) if pd.notna(row[4]) else None,
            "Xi": float(row[5]) if pd.notna(row[5]) else None,
            "Yi": float(row[6]) if pd.notna(row[6]) else None,
        })
    globals_data["circles"] = circles

    # === NON-CIRC SHAPES ===
    noncirc_df = xls.parse('non-circ')
    non_circ = list(noncirc_df.iloc[1:].dropna(subset=['Unnamed: 0']).apply(
        lambda row: {
            "X": float(row['Unnamed: 0']),
            "Y": float(row['Unnamed: 1']),
            "Movement": row['Unnamed: 2']
        }, axis=1))
    globals_data["non_circ"] = non_circ

    # === DISTRIBUTED LOADS ===
    # === DISTRIBUTED LOADS ===
    dload_df = xls.parse('dloads', header=None)
    dloads = []

    dload_data_blocks = [
        {"start_row": 3, "end_row": 13},   # Excel rows 4–13
        {"start_row": 16, "end_row": 26}   # Excel rows 17–26
    ]
    dload_block_starts = [1, 5, 9, 13]  # Corrected: Columns B:D, F:H, J:L, N:P

    for block in dload_data_blocks:
        for col in dload_block_starts:
            section = dload_df.iloc[block["start_row"]:block["end_row"], col:col+3].dropna()
            rows = section.apply(lambda row: {
                "X": float(row.iloc[0]),
                "Y": float(row.iloc[1]),
                "Normal": float(row.iloc[2])
            }, axis=1)
            block_points = list(rows)
            if block_points:
                dloads.append(block_points)
    globals_data["dloads"] = dloads

    # === REINFORCEMENT LINES ===
    reinforce_df = xls.parse('reinforce', header=None)
    reinforce_lines = []

    reinforce_data_blocks = [
        {"start_row": 3, "end_row": 13},   # Excel 4–13
        {"start_row": 16, "end_row": 26},  # Excel 17–26
        {"start_row": 29, "end_row": 39}   # Excel 30–39
    ]
    reinforce_block_starts = [1, 6, 11, 16]  # Columns B, G, L, Q

    for block in reinforce_data_blocks:
        for col in reinforce_block_starts:
            section = reinforce_df.iloc[block["start_row"]:block["end_row"], col:col+4].dropna()
            rows = section.apply(lambda row: {
                "X": float(row.iloc[0]),
                "Y": float(row.iloc[1]),
                "FL": float(row.iloc[2]),
                "FT": float(row.iloc[3])
            }, axis=1)
            line_points = list(rows)
            if line_points:
                reinforce_lines.append(line_points)
    globals_data["reinforce_lines"] = reinforce_lines

    # === STATIC GLOBALS ===
    globals_data["gamma_water"] = 62.4
    globals_data["tcrack_depth"] = 0.0
    globals_data["tcrack_water"] = 0.0

    return globals_data
