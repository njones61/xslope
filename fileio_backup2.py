import pandas as pd

def load_globals(filepath):
    """
    This function reads input data from various Excel sheets and parses it into
    structured components used throughout the slope stability analysis framework.
    It handles circular and non-circular failure surface data, reinforcement, piezometric
    lines, and distributed loads.

    Validation is enforced to ensure required geometry and material information is present:
    - Circular failure surface: must contain at least one valid row with Xo and Yo
    - Non-circular failure surface: required if no circular data is provided
    - Profile lines: must contain at least one valid set, and each line must have ≥ 2 points
    - Materials: must match the number of profile lines
    - Piezometric line: only included if it contains ≥ 2 valid rows
    - Distributed loads and reinforcement: each block must contain ≥ 2 valid entries

    Raises:
        ValueError: if required inputs are missing or inconsistent.

    Returns:
        dict: Parsed and validated global data structure for analysis
    """

    xls = pd.ExcelFile(filepath)
    globals_data = {}

    # === STATIC GLOBALS ===
    main_df = xls.parse('main', header=None)

    try:
        globals_data["gamma_water"] = float(main_df.iloc[15, 3])  # Excel row 16, column D
        globals_data["tcrack_depth"] = float(main_df.iloc[16, 3])  # Excel row 17, column D
        globals_data["tcrack_water"] = float(main_df.iloc[17, 3])  # Excel row 18, column D
        globals_data["k_seismic"] = float(main_df.iloc[18, 3])  # Excel row 19, column D
    except Exception as e:
        raise ValueError(f"Error reading static global values from 'main' tab: {e}")


    # === PROFILE LINES ===
    profile_df = xls.parse('profile', header=None)
    profile_lines = []

    profile_data_blocks = [
        {"header_row": 2, "data_start": 3, "data_end": 18},
        {"header_row": 20, "data_start": 21, "data_end": 36}
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
            data = profile_df.iloc[block["data_start"]:block["data_end"], [x_col, y_col]]
            data = data.dropna(how='all')
            if data.empty:
                continue
            if data.iloc[0].isna().any():
                continue
            coords = data.dropna().apply(lambda r: (float(r.iloc[0]), float(r.iloc[1])), axis=1).tolist()
            if len(coords) == 1:
                raise ValueError("Each profile line must contain at least two points.")
            if coords:
                profile_lines.append(coords)

    # === MATERIALS ===
    mat_df = xls.parse('mat', header=2)
    materials = []

    for _, row in mat_df.iterrows():
        try:
            gamma = float(row.iloc[1])  # Column B
            option = str(row.iloc[2]).strip().lower()  # Column C
            piezo = float(row.iloc[7]) if pd.notna(row.iloc[7]) else 1.0  # Column H

            if option == 'mc':
                if pd.isna(row.iloc[3]) or pd.isna(row.iloc[4]):  # D or E
                    continue  # missing c or phi
                c = float(row.iloc[3])
                phi = float(row.iloc[4])
                cp = 0
                r_elev = 0

            elif option == 'cp':
                if pd.isna(row.iloc[5]) or pd.isna(row.iloc[6]):  # F or G
                    continue  # missing cp or r_elev
                c = 0
                phi = 0
                cp = float(row.iloc[5])
                r_elev = float(row.iloc[6])

            else:
                continue  # unrecognized option

        except:
            continue  # skip if gamma, option, or piezo are invalid

        # Append parsed row
        materials.append({
            "gamma": gamma,
            "option": option,
            "c": c,
            "phi": phi,
            "cp": cp,
            "r_elev": r_elev,
            "piezo": piezo,
            "sigma_gamma": float(row.iloc[8]) if pd.notna(row.iloc[8]) else 0,
            "sigma_c": float(row.iloc[9]) if pd.notna(row.iloc[9]) else 0,
            "sigma_phi": float(row.iloc[10]) if pd.notna(row.iloc[10]) else 0,
            "sigma_cp": float(row.iloc[11]) if pd.notna(row.iloc[11]) else 0,
        })

    # === PIEZOMETRIC LINE ===
    piezo_df = xls.parse('piezo')
    piezo_data = piezo_df.iloc[2:].dropna(how='all')
    piezo_line = []
    if len(piezo_data.dropna()) >= 2:
        if len(piezo_data.dropna()) == 1:
            raise ValueError("Piezometric line must contain at least two points.")
        piezo_line = piezo_data.apply(lambda row: (float(row[0]), float(row[1])), axis=1).tolist()

    # === DISTRIBUTED LOADS ===
    dload_df = xls.parse('dloads', header=None)
    dloads = []
    dload_data_blocks = [
        {"start_row": 3, "end_row": 13},
        {"start_row": 16, "end_row": 26}
    ]
    dload_block_starts = [1, 5, 9, 13]
    for block in dload_data_blocks:
        for col in dload_block_starts:
            section = dload_df.iloc[block["start_row"]:block["end_row"], col:col+3].dropna(how='all')
            if len(section.dropna()) >= 2:
                rows = section.dropna().apply(lambda row: {
                    "X": float(row.iloc[0]),
                    "Y": float(row.iloc[1]),
                    "Normal": float(row.iloc[2])
                }, axis=1)
                block_points = list(rows)
                if len(block_points) == 1:
                    raise ValueError("Each distributed load must contain at least two points.")
                if block_points:
                    dloads.append(block_points)

    # === CIRCLES ===

    # Read the first 3 rows to get the max depth
    raw_df = xls.parse('circles', header=None)  # No header, get full sheet
    globals_data["max_depth"] = float(raw_df.iloc[0, 2])  # Excel C1 = row 0, column 2

    # Read the circles data starting from row 4 (index 3)
    circles_df = xls.parse('circles', header=3)
    raw = circles_df.dropna(subset=['Xo', 'Yo'], how='any')
    circles = []
    for _, row in raw.iterrows():
        circle = {
            "Xo": row['Xo'],
            "Yo": row['Yo'],
            "Option": row.get('Option', None),
            "Depth": row.get('Depth', None),
            "Xi": row.get('Xi', None),
            "Yi": row.get('Yi', None),
        }
        circles.append(circle)

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

    # === NON-CIRCULAR SURFACES ===
    noncirc_df = xls.parse('non-circ')
    non_circ = list(noncirc_df.iloc[1:].dropna(subset=['Unnamed: 0']).apply(
        lambda row: {
            "X": float(row['Unnamed: 0']),
            "Y": float(row['Unnamed: 1']),
            "Movement": row['Unnamed: 2']
        }, axis=1))

    # === REINFORCEMENT LINES ===
    reinforce_df = xls.parse('reinforce', header=None)
    reinforce_lines = []
    reinforce_data_blocks = [
        {"start_row": 3, "end_row": 13},
        {"start_row": 16, "end_row": 26},
        {"start_row": 29, "end_row": 39}
    ]
    reinforce_block_starts = [1, 6, 11, 16]
    for block in reinforce_data_blocks:
        for col in reinforce_block_starts:
            section = reinforce_df.iloc[block["start_row"]:block["end_row"], col:col+4].dropna(how='all')
            if len(section.dropna()) >= 2:
                rows = section.dropna().apply(lambda row: {
                    "X": float(row.iloc[0]),
                    "Y": float(row.iloc[1]),
                    "FL": float(row.iloc[2]),
                    "FT": float(row.iloc[3])
                }, axis=1)
                line_points = list(rows)
                if len(line_points) == 1:
                    raise ValueError("Each reinforcement line must contain at least two points.")
                if line_points:
                    reinforce_lines.append(line_points)

    # === VALIDATION ===
    circular = len(circles) > 0
    if not circular and len(globals_data.get('non_circ', [])) == 0:
        raise ValueError("Input must include either circular or non-circular surface data.")
    if not profile_lines:
        raise ValueError("Profile lines sheet is empty or invalid.")
    if not materials:
        raise ValueError("Materials sheet is empty.")
    if len(materials) != len(profile_lines):
        raise ValueError("Each profile line must have a corresponding material.")

    # Add everything to globals_data
    globals_data["profile_lines"] = profile_lines
    globals_data["materials"] = materials
    globals_data["piezo_line"] = piezo_line
    globals_data["circular"] = circular # True if circles are present
    globals_data["circles"] = circles
    globals_data["non_circ"] = non_circ
    globals_data["dloads"] = dloads
    globals_data["reinforce_lines"] = reinforce_lines

    return globals_data