import pandas as pd

def load_globals(filepath):
    """
    Parses an Excel file containing slopetools input and returns a structured dictionary.

    The Excel file is expected to contain several sheets with specific formats:
    - 'profile': Contains coordinate data for profile lines in multiple blocks.
    - 'mat': Contains material properties including strength parameters, pore pressure info, and uncertainty values.
    - 'piezo': Contains coordinates for the piezometric line.
    - 'circles': Contains slip circle definitions and maximum depth.
    - 'non-circ': Contains coordinates and movement options for non-circular failure surfaces.
    - 'dloads': Contains distributed load definitions in multiple blocks.
    - 'reinforce': Contains reinforcement line definitions in multiple blocks.

    Parameters:
        filepath (str): Path to the Excel file.

    Returns:
        dict: A dictionary containing the parsed global data with the following keys:
            - "profile_lines": List of profile line segments (list of (x, y) tuples).
            - "materials": List of material property dictionaries.
            - "piezo_line": List of (x, y) tuples representing the piezometric line.
            - "max_depth": Maximum analysis depth from the 'circles' sheet.
            - "circles": List of dictionaries defining circular slip surfaces.
            - "non_circ": List of dictionaries defining non-circular surface points.
            - "dloads": List of distributed load blocks, each a list of (x, y, normal) dicts.
            - "reinforce_lines": List of reinforcement lines, each a list of (x, y, FL, FT) dicts.
            - "gamma_water": Static water unit weight (default: 62.4).
            - "tcrack_depth": Top crack depth (default: 0.0).
            - "tcrack_water": Top crack water depth (default: 0.0).

    Notes:
        - The function assumes fixed layouts for rows and columns in each sheet.
        - Malformed or incomplete rows are skipped.
        - Sheet names must match exactly as described above.

    Raises:
        KeyError: If required columns or sheets are missing.
        ValueError: If critical data cannot be parsed into the expected format.
    """

    xls = pd.ExcelFile(filepath)
    globals_data = {}

    # === STATIC GLOBALS ===
    main_df = xls.parse('main', header=None)

    globals_data["gamma_water"] = float(main_df.iloc[15, 3]) # Excel row 16, column D
    globals_data["tcrack_depth"] = float(main_df.iloc[16, 3]) # Excel row 17, column D
    globals_data["tcrack_water"] = float(main_df.iloc[17, 3]) # Excel row 18, column D
    globals_data["k_seismic"] = float(main_df.iloc[18, 3])    # Excel row 19, column D

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
    mat_df = xls.parse('mat', header=2)

    materials = []
    for _, row in mat_df.iterrows():
        if row.iloc[1:12].dropna().empty:
            continue  # Skip row if columns B–L are all empty
        try:
            gamma = float(row.iloc[1])
            option = row.iloc[2]
            c = float(row.iloc[3])
            phi = float(row.iloc[4])
            cp = float(row.iloc[5])
            r_elev = float(row.iloc[6])
        except:
            continue  # Skip rows with missing or invalid data
        materials.append({
            "gamma": gamma,
            "option": option,
            "c": c,
            "phi": phi,
            "cp": cp,
            "r_elev": r_elev,
            "piezo": float(row['Piezo']),
            "sigma_gamma": float(row.iloc[7]) if pd.notna(row.iloc[7]) else 0,
            "sigma_c": float(row.iloc[8]) if pd.notna(row.iloc[8]) else 0,
            "sigma_phi": float(row.iloc[9]) if pd.notna(row.iloc[9]) else 0,
            "sigma_cp": float(row.iloc[10]) if pd.notna(row.iloc[10]) else 0,
        })

    globals_data["materials"] = materials

    # === PIEZOMETRIC LINE ===
    piezo_df = xls.parse('piezo')
    piezo_line = list(piezo_df.iloc[2:].dropna().apply(lambda row: (float(row.iloc[0]), float(row.iloc[1])), axis=1))
    globals_data["piezo_line"] = piezo_line

    # === CIRCLES ===
    circles_df = xls.parse('circles')
    globals_data["max_depth"] = float(circles_df.iloc[0, 2])
    circles_rows = circles_df.iloc[3:]
    circles_rows = circles_rows.dropna(subset=circles_df.columns[[1, 2]])
    circles = []
    for _, row in circles_rows.iterrows():
        circles.append({
            "Xo": float(row.iloc[1]),
            "Yo": float(row.iloc[2]),
            "Option": row.iloc[3],
            "Depth": float(row.iloc[4]) if pd.notna(row.iloc[4]) else None,
            "Xi": float(row.iloc[5]) if pd.notna(row.iloc[5]) else None,
            "Yi": float(row.iloc[6]) if pd.notna(row.iloc[6]) else None,
        })
    globals_data["circles"] = circles

    # For each circle, fill in the radius and depth values depending on the circle option
    for circle in circles:
        if circle["Option"] == "Depth":
            yo = circle["Yo"]
            depth = circle["Depth"]
            circle["R"] = yo - depth
        elif circle["Option"] == "Intercept":
            xi = circle["Xi"]
            yi = circle["Yi"]
            xo = circle["Xo"]
            yo = circle["Yo"]
            r = ((xi - xo) ** 2 + (yi - yo) ** 2) ** 0.5
            circle["R"] = r
            circle["Depth"] = yo - r
        elif circle["Option"] == "Radius":
            yo = circle["Yo"]
            r = circle["R"]
            circle["Depth"] = yo - r
        else:
            raise ValueError(f"Unknown option '{circle['Option']}' for circles.")

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

    # === CHECK FOR CIRCULAR OR NON-CIRCULAR ===
    circular = True
    if len(circles) == 0:
        circular = False
        if len(non_circ) == 0:
            raise ValueError("Input must include either circular or non-circular surface data.")
    globals_data["circular"] = circular
    print(f"Input is {'circular' if circular else 'non-circular'}.")

    # === CHECK FOR EMPTY PROFILE LINES ===
    if not profile_lines:
        raise ValueError("Profile lines sheet is empty.")

    # === CHECK FOR EMPTY MATERIAL TABLE OR COUNT MISMATCH ===
    if not materials:
        raise ValueError("Materials sheet is empty.")
    if len(materials) != len(profile_lines):
        raise ValueError("Number of materials does not match number of profile lines.")

    return globals_data
