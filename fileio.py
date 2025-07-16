import pandas as pd
from shapely.geometry import LineString, Point
import numpy as np
import pickle
import os

# Import functions for seep analysis
from mesh import import_mesh_from_json

def build_ground_surface(profile_lines):
    """
    Constructs the topmost ground surface LineString from a set of profile lines.

    The function finds the highest elevation at each x-coordinate across all profile lines,
    which represents the true ground surface.

    Parameters:
        profile_lines (list of list of tuple): A list of profile lines, each represented
            as a list of (x, y) coordinate tuples.

    Returns:
        shapely.geometry.LineString: A LineString of the top surface, or an empty LineString
        if fewer than two valid points are found.
    """
    
    if not profile_lines:
        return LineString([])
    
    # Step 1: Gather all points from all profile lines
    all_points = []
    for line in profile_lines:
        all_points.extend(line)
    
    # Step 2: Group points by x-coordinate and find the highest y for each x
    x_groups = {}
    for x, y in all_points:
        if x not in x_groups:
            x_groups[x] = y
        else:
            x_groups[x] = max(x_groups[x], y)
    
    # Step 3: For each candidate point, check if any profile line is above it
    ground_surface_points = []
    for x, y in sorted(x_groups.items()):
        # Create a vertical line at this x-coordinate
        vertical_line = LineString([(x, y - 1000), (x, y + 1000)])
        
        # Check intersections with all profile lines
        is_topmost = True
        for profile_line in profile_lines:
            line = LineString(profile_line)
            if line.length == 0:
                continue
            
            # Find intersection with this profile line
            intersection = line.intersection(vertical_line)
            if not intersection.is_empty:
                # Get the y-coordinate of the intersection
                if hasattr(intersection, 'y'):
                    # Single point intersection
                    if intersection.y > y + 1e-6:  # Allow small numerical tolerance
                        is_topmost = False
                        break
                elif hasattr(intersection, 'geoms'):
                    # Multiple points or line intersection
                    for geom in intersection.geoms:
                        if hasattr(geom, 'y') and geom.y > y + 1e-6:
                            is_topmost = False
                            break
                    if not is_topmost:
                        break
        
        if is_topmost:
            ground_surface_points.append((x, y))
    
    # Ensure we have at least 2 points
    if len(ground_surface_points) < 2:
        return LineString([])
    
    return LineString(ground_surface_points)



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
        template_version = main_df.iloc[4, 3]  # Excel row 5, column D
        gamma_water = float(main_df.iloc[18, 3])  # Excel row 19, column D
        tcrack_depth = float(main_df.iloc[19, 3])  # Excel row 20, column D
        tcrack_water = float(main_df.iloc[20, 3])  # Excel row 21, column D
        k_seismic = float(main_df.iloc[21, 3])  # Excel row 22, column D
    except Exception as e:
        raise ValueError(f"Error reading static global values from 'main' tab: {e}")


    # === PROFILE LINES ===
    profile_df = xls.parse('profile', header=None)

    max_depth = float(profile_df.iloc[1, 1])  # Excel B2 = row 1, column 1

    profile_lines = []

    profile_data_blocks = [
        {"header_row": 4, "data_start": 5, "data_end": 20},
        {"header_row": 22, "data_start": 23, "data_end": 38}
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

    # === BUILD GROUND SURFACE FROM PROFILE LINES ===

    ground_surface = build_ground_surface(profile_lines)

    # === BUILD TENSILE CRACK LINE ===

    tcrack_surface = None
    if tcrack_depth > 0:
        tcrack_surface = LineString([(x, y - tcrack_depth) for (x, y) in ground_surface.coords])

    # === MATERIALS (Optimized Parsing) ===
    mat_df = xls.parse('mat', header=2)  # header=2 because the header row is row 3 in Excel
    materials = []

    # Only process rows 4-15 (Excel), which are 0-indexed 0-11 in pandas 
    for _, row in mat_df.iloc[0:12].iterrows():
        # Check if the row is blank (columns 2-22, which are indices 1-21)
        if row.iloc[1:22].isna().all():
            continue
        materials.append({
            "name": row.get('name', ''),
            "gamma": float(row.get('g', 0) or 0),
            "option": str(row.get('option', '')).strip().lower(),
            "c": float(row.get('c', 0) or 0),
            "phi": float(row.get('f', 0) or 0),
            "cp": float(row.get('cp', 0) or 0),
            "r_elev": float(row.get('r-elev', 0) or 0),
            "d": float(row.get('d', 0)) if pd.notna(row.get('d')) else 0,
            "psi": float(row.get('ψ', 0)) if pd.notna(row.get('ψ')) else 0,
            "u": str(row.get('u', 'none')).strip().lower(),
            "sigma_gamma": float(row.get('s(g)', 0) or 0),
            "sigma_c": float(row.get('s(c)', 0) or 0),
            "sigma_phi": float(row.get('s(f)', 0) or 0),
            "sigma_cp": float(row.get('s(cp)', 0) or 0),
            "sigma_d": float(row.get('s(d)', 0) or 0),
            "sigma_psi": float(row.get('s(ψ)', 0) or 0),
            "k1": float(row.get('k1', 0) or 0),
            "k2": float(row.get('k2', 0) or 0),
            "alpha": float(row.get('alpha', 0) or 0),
            "kr0" : float(row.get('kr0', 0) or 0),
            "h0" : float(row.get('h0', 0) or 0)
        })

    # === SEEPAGE ANALYSIS FILES ===
    # Check if any materials use seepage analysis for pore pressure
    has_seep_materials = any(material["u"] == "seep" for material in materials)
    
    seep_mesh = None
    seep_u = None
    seep_u2 = None
    
    if has_seep_materials:
        # Read seepage file names directly from Excel cells L19, L20, L21
        try:
            # Read the 'mat' sheet directly without header parsing
            mat_raw_df = xls.parse('mat', header=None)
            
            # L19 = row 18, column 11 (0-indexed)
            mesh_filename = str(mat_raw_df.iloc[18, 11]).strip()  # L19
            solution1_filename = str(mat_raw_df.iloc[19, 11]).strip()  # L20  
            solution2_filename = str(mat_raw_df.iloc[20, 11]).strip()  # L21
            
            # Validate required files
            if not mesh_filename or mesh_filename.lower() == 'nan':
                raise ValueError("CRITICAL ERROR: Mesh filename is required when using 'seep' pore pressure option but is blank in cell L19.")
            if not solution1_filename or solution1_filename.lower() == 'nan':
                raise ValueError("CRITICAL ERROR: Solution1 filename is required when using 'seep' pore pressure option but is blank in cell L20.")
            
            # Load mesh file
            if not os.path.exists(mesh_filename):
                raise ValueError(f"CRITICAL ERROR: Mesh file '{mesh_filename}' not found.")
            seep_mesh = import_mesh_from_json(mesh_filename)
            
            # Load solution1 file
            if not os.path.exists(solution1_filename):
                raise ValueError(f"CRITICAL ERROR: Solution1 file '{solution1_filename}' not found.")
            solution1_df = pd.read_csv(solution1_filename)
            # Skip the last row which contains the total flowrate comment
            solution1_df = solution1_df.iloc[:-1]
            seep_u = solution1_df["u"].to_numpy()
            
            # Load solution2 file if provided
            if solution2_filename and solution2_filename.lower() != 'nan':
                if not os.path.exists(solution2_filename):
                    raise ValueError(f"CRITICAL ERROR: Solution2 file '{solution2_filename}' not found.")
                solution2_df = pd.read_csv(solution2_filename)
                # Skip the last row which contains the total flowrate comment
                solution2_df = solution2_df.iloc[:-1]
                seep_u2 = solution2_df["u"].to_numpy()
                
        except Exception as e:
            if "CRITICAL ERROR" in str(e):
                raise e
            else:
                raise ValueError(f"Error reading seepage files: {e}")

    # === PIEZOMETRIC LINE ===
    piezo_df = xls.parse('piezo')
    piezo_line = []
    piezo_line2 = []

    # Read all data once (rows 4-18)
    piezo_data = piezo_df.iloc[2:18].dropna(how='all')
    
    if len(piezo_data) >= 2:
        # Extract first table (A4:B18) - columns 0 and 1
        try:
            piezo_data1 = piezo_data.dropna(subset=[piezo_data.columns[0], piezo_data.columns[1]], how='all')
            if len(piezo_data1) < 2:
                raise ValueError("First piezometric line must contain at least two points.")
            piezo_line = piezo_data1.apply(lambda row: (float(row.iloc[0]), float(row.iloc[1])), axis=1).tolist()
        except Exception:
            raise ValueError("Invalid first piezometric line format.")

        # Extract second table (D4:E18) - columns 3 and 4
        try:
            piezo_data2 = piezo_data.dropna(subset=[piezo_data.columns[3], piezo_data.columns[4]], how='all')
            if len(piezo_data2) < 2:
                raise ValueError("Second piezometric line must contain at least two points.")
            piezo_line2 = piezo_data2.apply(lambda row: (float(row.iloc[3]), float(row.iloc[4])), axis=1).tolist()
        except Exception:
            # If second table reading fails, just leave piezo_line2 as empty list
            piezo_line2 = []
    elif len(piezo_data) == 1:
        raise ValueError("Piezometric line must contain at least two points.")

    # === DISTRIBUTED LOADS ===
    dload_df = xls.parse('dloads', header=None)
    dloads = []
    dloads2 = []
    dload_data_blocks = [
        {"start_row": 3, "end_row": 13},
        {"start_row": 16, "end_row": 26}
    ]
    dload_block_starts = [1, 5, 9, 13]

    for block_idx, block in enumerate(dload_data_blocks):
        for col in dload_block_starts:
            section = dload_df.iloc[block["start_row"]:block["end_row"], col:col + 3]
            section = section.dropna(how='all')
            section = section.dropna(subset=[col, col + 1], how='any')
            if len(section) >= 2:
                try:
                    block_points = section.apply(
                        lambda row: {
                            "X": float(row.iloc[0]),
                            "Y": float(row.iloc[1]),
                            "Normal": float(row.iloc[2])
                        }, axis=1).tolist()
                    if block_idx == 0:
                        dloads.append(block_points)
                    else:
                        dloads2.append(block_points)
                except:
                    raise ValueError("Invalid data format in distributed load block.")
            elif len(section) == 1:
                raise ValueError("Each distributed load block must contain at least two points.")

    # === CIRCLES ===

    # Read the first 3 rows to get the max depth
    raw_df = xls.parse('circles', header=None)  # No header, get full sheet

    # Read the circles data starting from row 2 (index 1)
    circles_df = xls.parse('circles', header=1)
    raw = circles_df.dropna(subset=['Xo', 'Yo'], how='any')
    circles = []
    for _, row in raw.iterrows():
        Xo = row['Xo']
        Yo = row['Yo']
        Option = row.get('Option', None)
        Depth = row.get('Depth', None)
        Xi = row.get('Xi', None)
        Yi = row.get('Yi', None)
        R = row.get('R', None)
        # For each circle, fill in the radius and depth values depending on the circle option
        if Option == 'Depth':
            R = Yo - Depth
        elif Option == 'Intercept':
            R = ((Xi - Xo) ** 2 + (Yi - Yo) ** 2) ** 0.5
            Depth = Yo - R
        elif Option == 'Radius':
            Depth = Yo - R
        else:
            raise ValueError(f"Unknown option '{Option}' for circles.")
        circle = {
            "Xo": Xo,
            "Yo": Yo,
            "Depth": Depth,
            "R": R,
        }
        circles.append(circle)

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
            section = reinforce_df.iloc[block["start_row"]:block["end_row"], col:col + 4]
            section = section.dropna(how='all')
            section = section.dropna(subset=[col, col + 1], how='any')
            if len(section) >= 2:
                try:
                    line_points = section.apply(
                        lambda row: {
                            "X": float(row.iloc[0]),
                            "Y": float(row.iloc[1]),
                            "FL": float(row.iloc[2]),
                            "FT": float(row.iloc[3])
                        }, axis=1).tolist()
                    reinforce_lines.append(line_points)
                except:
                    raise ValueError("Invalid data format in reinforcement block.")
            elif len(section) == 1:
                raise ValueError("Each reinforcement line must contain at least two points.")


    # === SEEPAGE ANALYSIS BOUNDARY CONDITIONS ===
    seep_df = xls.parse('seep bc', header=None)
    seepage_bc = {"specified_heads": [], "exit_face": []}

    # Specified Head #1
    head1 = seep_df.iloc[2, 2] if seep_df.shape[1] > 2 and seep_df.shape[0] > 2 else None
    coords1 = []
    for i in range(4, 12):  # rows 5-12 (0-indexed 4-11)
        if i >= seep_df.shape[0]:
            break
        x = seep_df.iloc[i, 1] if seep_df.shape[1] > 1 else None
        y = seep_df.iloc[i, 2] if seep_df.shape[1] > 2 else None
        if pd.notna(x) and pd.notna(y):
            coords1.append((float(x), float(y)))
    if head1 is not None and coords1:
        seepage_bc["specified_heads"].append({"head": float(head1), "coords": coords1})

    # Specified Head #2
    head2 = seep_df.iloc[2, 5] if seep_df.shape[1] > 5 and seep_df.shape[0] > 2 else None
    coords2 = []
    for i in range(4, 12):
        if i >= seep_df.shape[0]:
            break
        x = seep_df.iloc[i, 4] if seep_df.shape[1] > 4 else None
        y = seep_df.iloc[i, 5] if seep_df.shape[1] > 5 else None
        if pd.notna(x) and pd.notna(y):
            coords2.append((float(x), float(y)))
    if head2 is not None and coords2:
        seepage_bc["specified_heads"].append({"head": float(head2), "coords": coords2})

    # Specified Head #3
    head3 = seep_df.iloc[2, 8] if seep_df.shape[1] > 8 and seep_df.shape[0] > 2 else None
    coords3 = []
    for i in range(4, 12):
        if i >= seep_df.shape[0]:
            break
        x = seep_df.iloc[i, 7] if seep_df.shape[1] > 7 else None
        y = seep_df.iloc[i, 8] if seep_df.shape[1] > 8 else None
        if pd.notna(x) and pd.notna(y):
            coords3.append((float(x), float(y)))
    if head3 is not None and coords3:
        seepage_bc["specified_heads"].append({"head": float(head3), "coords": coords3})

    # Exit Face
    exit_coords = []
    for i in range(15, 23):  # rows 16-23 (0-indexed 15-22)
        if i >= seep_df.shape[0]:
            break
        x = seep_df.iloc[i, 1] if seep_df.shape[1] > 1 else None
        y = seep_df.iloc[i, 2] if seep_df.shape[1] > 2 else None
        if pd.notna(x) and pd.notna(y):
            exit_coords.append((float(x), float(y)))
    seepage_bc["exit_face"] = exit_coords

    # === VALIDATION ===
 
    circular = len(circles) > 0
    if not circular and len(non_circ) == 0:
        raise ValueError("Input must include either circular or non-circular surface data.")
    if not profile_lines:
        raise ValueError("Profile lines sheet is empty or invalid.")
    if not materials:
        raise ValueError("Materials sheet is empty.")
    if len(materials) != len(profile_lines):
        raise ValueError("Each profile line must have a corresponding material. You have " + str(len(materials)) + " materials and " + str(len(profile_lines)) + " profile lines.")
        

    # Add everything to globals_data
    globals_data["template_version"] = template_version
    globals_data["gamma_water"] = gamma_water
    globals_data["tcrack_depth"] = tcrack_depth
    globals_data["tcrack_water"] = tcrack_water
    globals_data["k_seismic"] = k_seismic
    globals_data["max_depth"] = max_depth
    globals_data["profile_lines"] = profile_lines
    globals_data["ground_surface"] = ground_surface
    globals_data["tcrack_surface"] = tcrack_surface
    globals_data["materials"] = materials
    globals_data["piezo_line"] = piezo_line
    globals_data["piezo_line2"] = piezo_line2
    globals_data["circular"] = circular # True if circles are present
    globals_data["circles"] = circles
    globals_data["non_circ"] = non_circ
    globals_data["dloads"] = dloads
    globals_data["dloads2"] = dloads2
    globals_data["reinforce_lines"] = reinforce_lines
    globals_data["seepage_bc"] = seepage_bc
    
    # Add seepage data if available
    if has_seep_materials:
        globals_data["seep_mesh"] = seep_mesh
        globals_data["seep_u"] = seep_u
        if seep_u2 is not None:
            globals_data["seep_u2"] = seep_u2

    return globals_data

def save_data_to_pickle(data, filepath):
    """
    Save a data object to a pickle file.
    
    This function serializes the data object and saves it to the specified filepath.
    Useful for saving processed data from Excel templates for later use.
    
    Parameters:
        data: The data object to save (typically a dictionary from load_globals)
        filepath (str): The file path where the pickle file should be saved
        
    Returns:
        None
        
    Raises:
        IOError: If the file cannot be written
        PickleError: If the data cannot be serialized
    """
    try:
        with open(filepath, 'wb') as f:
            pickle.dump(data, f)
    except Exception as e:
        raise IOError(f"Failed to save data to pickle file '{filepath}': {e}")


def load_data_from_pickle(filepath):
    """
    Load a data object from a pickle file.
    
    This function deserializes a data object from the specified pickle file.
    Useful for loading previously saved data without re-processing Excel templates.
    
    Parameters:
        filepath (str): The file path of the pickle file to load
        
    Returns:
        The deserialized data object (typically a dictionary)
        
    Raises:
        FileNotFoundError: If the pickle file doesn't exist
        IOError: If the file cannot be read
        PickleError: If the data cannot be deserialized
    """
    try:
        with open(filepath, 'rb') as f:
            data = pickle.load(f)
        return data
    except FileNotFoundError:
        raise FileNotFoundError(f"Pickle file not found: '{filepath}'")
    except Exception as e:
        raise IOError(f"Failed to load data from pickle file '{filepath}': {e}")