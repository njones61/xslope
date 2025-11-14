#### PLACE TO PUT MISC DEBUGGING CODE ####


def compare_python_fortran_nodal_results(python_head, python_q, fortran_out_path, node_offset=1, verbose=True, nbc=None):
    """
    Compare Python and FORTRAN nodal heads and flows and save results to Excel.
    - python_head: numpy array of nodal heads (Python, index 0 = node 1 in FORTRAN)
    - python_q: numpy array of nodal flows (same indexing)
    - fortran_out_path: path to FORTRAN .out file
    - node_offset: 1 if FORTRAN nodes are 1-based, 0 if 0-based
    - verbose: if True, print summary statistics
    - nbc: array of boundary condition types (0=interior, 1=fixed head, 2=exit face)
    """
    import re
    import numpy as np
    import pandas as pd
    from pathlib import Path

    # Debug prints for nbc array
    print("\n=== Debug: Boundary Conditions ===")
    print(f"nbc is None: {nbc is None}")
    if nbc is not None:
        print(f"nbc type: {type(nbc)}")
        print(f"nbc shape: {nbc.shape if hasattr(nbc, 'shape') else len(nbc)}")
        print(f"nbc dtype: {nbc.dtype if hasattr(nbc, 'dtype') else type(nbc[0])}")
        print(f"nbc unique values: {np.unique(nbc)}")
        print(f"nbc first 10 values: {nbc[:10]}")
        print(f"Number of interior nodes (0): {np.sum(nbc == 0)}")
        print(f"Number of fixed head nodes (1): {np.sum(nbc == 1)}")
        print(f"Number of exit face nodes (2): {np.sum(nbc == 2)}")

    # Read the FORTRAN output file and find the relevant section
    with open(fortran_out_path, 'r') as f:
        lines = f.readlines()

    # Find the start of the 'Nodal Flows and Heads' section
    start_idx = None
    for i, line in enumerate(lines):
        if 'Nodal Flows and Heads' in line:
            start_idx = i
            break
    if start_idx is None:
        raise ValueError('Could not find "Nodal Flows and Heads" section in FORTRAN output.')

    # Skip header lines to the data
    data_start = start_idx + 5  # 5 lines after header is usually where data starts
    nodal_data = []
    for line in lines[data_start:]:
        if not line.strip():
            continue
        # Stop if we hit a non-data line (e.g., next section)
        if re.match(r'\s*\d+\s+\d', line) is None and re.match(r'\s*\d+\s+\d', line.strip()) is None:
            # If line doesn't start with a node number, break
            if not re.match(r'\s*\d+', line):
                break
        # Parse node, head, percent, (optional) flow
        m = re.match(r'\s*(\d+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s*%?\s*([\d.Ee+-]*)', line)
        if m:
            node = int(m.group(1))
            head = float(m.group(2))
            # percent = float(m.group(3))  # Not used
            flow_str = m.group(4)
            flow = float(flow_str) if flow_str.strip() else np.nan
            nodal_data.append((node, head, flow))
        else:
            # If the line doesn't match, stop parsing
            break

    # Convert to arrays
    fortran_nodes = np.array([n[0] for n in nodal_data], dtype=int)
    fortran_head = np.array([n[1] for n in nodal_data], dtype=float)
    fortran_flow = np.array([n[2] for n in nodal_data], dtype=float)

    # Align indices: assume node 1 in FORTRAN is index 0 in Python
    n_compare = min(len(fortran_head), len(python_head))
    head_diff = python_head[:n_compare] - fortran_head[:n_compare]
    q_diff = python_q[:n_compare] - fortran_flow[:n_compare]

    # Create BC type mapping
    bc_type_map = {0: 'Interior', 1: 'Fixed Head', 2: 'Exit Face'}
    bc_types = np.array(['Unknown'] * n_compare)
    if nbc is not None:
        print(f"\nDebug: Creating BC types array")
        print(f"n_compare: {n_compare}")
        print(f"nbc length: {len(nbc)}")
        # Convert nbc to numpy array if it isn't already
        nbc = np.asarray(nbc)
        # Create BC types array
        bc_types = np.array([bc_type_map.get(int(nbc[i]), 'Unknown') for i in range(n_compare)])
        print(f"First 10 BC types: {bc_types[:10]}")
        print(f"Unique BC types: {np.unique(bc_types)}")

    # Create DataFrame for comparison
    df = pd.DataFrame({
        'Node': np.arange(node_offset, n_compare + node_offset),
        'BC_Type': bc_types,
        'Python_Head': python_head[:n_compare],
        'FORTRAN_Head': fortran_head[:n_compare],
        'Head_Diff': head_diff,
        'Abs(Head_Diff)': np.abs(head_diff),
        'Python_Q': python_q[:n_compare],
        'FORTRAN_Q': fortran_flow[:n_compare],
        'Q_Diff': q_diff,
        'Abs(Q_Diff)': np.abs(q_diff)
    })

    # Add summary statistics
    summary_stats = {
        'Metric': [
            'Max abs(head diff)',
            'Mean abs(head diff)',
            'Std abs(head diff)',
            'Max abs(q diff)',
            'Mean abs(q diff)',
            'Std abs(q diff)'
        ],
        'Value': [
            np.max(np.abs(head_diff)),
            np.mean(np.abs(head_diff)),
            np.std(head_diff),
            np.max(np.abs(q_diff[~np.isnan(fortran_flow[:n_compare])])),
            np.mean(np.abs(q_diff[~np.isnan(fortran_flow[:n_compare])])),
            np.std(q_diff[~np.isnan(fortran_flow[:n_compare])])
        ]
    }
    summary_df = pd.DataFrame(summary_stats)

    # Add BC-specific statistics if nbc is provided
    if nbc is not None:
        bc_stats = []
        for bc_type in [0, 1, 2]:
            mask = nbc[:n_compare] == bc_type
            if np.any(mask):
                bc_name = bc_type_map[bc_type]
                # Head statistics
                head_mask = mask
                if np.any(head_mask):
                    bc_stats.extend([
                        [f'{bc_name} - Max abs(head diff)', np.max(np.abs(head_diff[head_mask]))],
                        [f'{bc_name} - Mean abs(head diff)', np.mean(np.abs(head_diff[head_mask]))]
                    ])
                
                # Flow statistics (only for nodes with valid FORTRAN flow data)
                flow_mask = mask & ~np.isnan(fortran_flow[:n_compare])
                if np.any(flow_mask):
                    bc_stats.extend([
                        [f'{bc_name} - Max abs(q diff)', np.max(np.abs(q_diff[flow_mask]))],
                        [f'{bc_name} - Mean abs(q diff)', np.mean(np.abs(q_diff[flow_mask]))]
                    ])
        
        if bc_stats:
            bc_summary_df = pd.DataFrame(bc_stats, columns=['Metric', 'Value'])
            summary_df = pd.concat([summary_df, bc_summary_df], ignore_index=True)

    # Create Excel writer
    output_path = Path(fortran_out_path).with_suffix('.xlsx')
    with pd.ExcelWriter(output_path) as writer:
        # Write main comparison data
        df.to_excel(writer, sheet_name='Nodal Comparison', index=False)
        
        # Write summary statistics
        summary_df.to_excel(writer, sheet_name='Summary Statistics', index=False)
        
        # Auto-adjust column widths
        for sheet_name in writer.sheets:
            worksheet = writer.sheets[sheet_name]
            for idx, col in enumerate(df.columns):
                max_length = max(
                    df[col].astype(str).apply(len).max(),
                    len(col)
                )
                worksheet.column_dimensions[chr(65 + idx)].width = max_length + 2

    if verbose:
        print(f"\nComparison results saved to: {output_path}")
        print("\n=== Python vs FORTRAN Nodal Head Comparison ===")
        print(f"Max abs(head diff): {np.max(np.abs(head_diff)):.4e}")
        print(f"Mean abs(head diff): {np.mean(np.abs(head_diff)):.4e}")
        print(f"Std abs(head diff): {np.std(head_diff):.4e}")

        print("\n=== Python vs FORTRAN Nodal Flow Comparison ===")
        # Only compare where FORTRAN flow is not nan
        valid = ~np.isnan(fortran_flow[:n_compare])
        if np.any(valid):
            print(f"Max abs(q diff): {np.max(np.abs(q_diff[valid])):.4e}")
            print(f"Mean abs(q diff): {np.mean(np.abs(q_diff[valid])):.4e}")
            print(f"Std abs(q diff): {np.std(q_diff[valid]):.4e}")
        else:
            print("No FORTRAN flow data to compare.")

        if nbc is not None:
            print("\n=== Boundary Condition Statistics ===")
            for bc_type in [0, 1, 2]:
                mask = nbc[:n_compare] == bc_type
                if np.any(mask):
                    bc_name = bc_type_map[bc_type]
                    print(f"\n{bc_name} nodes:")
                    print(f"  Count: {np.sum(mask)}")
                    # Head statistics
                    print(f"  Max abs(head diff): {np.max(np.abs(head_diff[mask])):.4e}")
                    print(f"  Mean abs(head diff): {np.mean(np.abs(head_diff[mask])):.4e}")
                    # Flow statistics (only for nodes with valid FORTRAN flow data)
                    valid_q = mask & ~np.isnan(fortran_flow[:n_compare])
                    if np.any(valid_q):
                        print(f"  Max abs(q diff): {np.max(np.abs(q_diff[valid_q])):.4e}")
                        print(f"  Mean abs(q diff): {np.mean(np.abs(q_diff[valid_q])):.4e}")
                    else:
                        print("  No valid flow data for comparison")

    return head_diff, q_diff

def plot_kr_field(coords, elements, kr_vals, title='Kr Field'):
    """
    Plot the kr field at element centroids.
    """
    import numpy as np
    centroids = np.mean(coords[elements], axis=1)
    plt.figure(figsize=(10, 4))
    sc = plt.scatter(centroids[:, 0], centroids[:, 1], c=kr_vals, cmap='viridis', s=30, edgecolor='k')
    plt.colorbar(sc, label='k_r')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(title)
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    plt.show()

def plot_nodal_flows(coords, q, title='Nodal Flow Values', show_values=False, cmap='RdBu_r'):
    """
    Plot nodal flow values for debugging.
    
    Parameters:
    -----------
    coords : ndarray
        Node coordinates array (n_nodes x 2)
    q : ndarray
        Nodal flow values array (n_nodes)
    title : str, optional
        Plot title
    show_values : bool, optional
        Whether to show the actual flow values on the plot
    cmap : str, optional
        Colormap to use for the scatter plot
    """
    plt.figure(figsize=(10, 8))
    
    # Create scatter plot of nodes colored by flow value
    scatter = plt.scatter(coords[:, 0], coords[:, 1], 
                         c=q, cmap=cmap, 
                         s=100, edgecolor='k')
    
    # Add colorbar
    cbar = plt.colorbar(scatter)
    cbar.set_label('Flow Value')
    
    # Add node numbers and flow values if requested
    if show_values:
        for i, (x, y) in enumerate(coords):
            plt.text(x, y, f'Node {i+1}\n{q[i]:.2e}', 
                    ha='center', va='center',
                    fontsize=8)
    
    plt.title(title)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid(True)
    plt.axis('equal')
    
    return plt.gcf()

def classify_nodes_by_phreatic_surface(coords, head):
    """
    Returns a boolean mask: True if node is above the phreatic surface (unsaturated), False otherwise.
    """
    return head < coords[:, 1]

def debug_nodal_flows_above_phreatic(coords, head, q, title='Nodal Flows vs Phreatic Surface'):
    above_mask = classify_nodes_by_phreatic_surface(coords, head)
    below_mask = ~above_mask

    # Print summary statistics
    print(f'Nodes above phreatic surface: {np.sum(above_mask)}')
    print(f'Nodes below phreatic surface: {np.sum(below_mask)}')
    print('--- Above phreatic surface ---')
    print(f'Max |q|: {np.max(np.abs(q[above_mask])):.3e}')
    print(f'Mean |q|: {np.mean(np.abs(q[above_mask])):.3e}')
    print(f'Min |q|: {np.min(np.abs(q[above_mask])):.3e}')
    print('--- Below phreatic surface ---')
    print(f'Max |q|: {np.max(np.abs(q[below_mask])):.3e}')
    print(f'Mean |q|: {np.mean(np.abs(q[below_mask])):.3e}')
    print(f'Min |q|: {np.min(np.abs(q[below_mask])):.3e}')

    # Plot
    plt.figure(figsize=(10, 8))
    plt.scatter(coords[below_mask, 0], coords[below_mask, 1], c=q[below_mask], cmap='Blues', label='Below Phreatic', s=60)
    plt.scatter(coords[above_mask, 0], coords[above_mask, 1], c=q[above_mask], cmap='Reds', label='Above Phreatic', s=60, marker='^')
    plt.colorbar(label='Nodal q')
    plt.legend()
    plt.title(title)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid(True)
    plt.axis('equal')
    plt.show()

def debug_kr_above_phreatic(coords, head, kr, title='kr vs Phreatic Surface'):
    """
    Plot and summarize kr values above and below the phreatic surface.
    coords: (n_nodes, 2) array of node coordinates
    head: (n_nodes,) array of nodal heads
    kr: (n_nodes,) array of nodal kr values
    """
    import numpy as np
    import matplotlib.pyplot as plt
    above_mask = classify_nodes_by_phreatic_surface(coords, head)
    below_mask = ~above_mask

    # Print summary statistics
    print(f'Nodes above phreatic surface: {np.sum(above_mask)}')
    print(f'Nodes below phreatic surface: {np.sum(below_mask)}')
    print('--- Above phreatic surface ---')
    print(f'Max kr: {np.max(kr[above_mask]):.3e}')
    print(f'Mean kr: {np.mean(kr[above_mask]):.3e}')
    print(f'Min kr: {np.min(kr[above_mask]):.3e}')
    print('--- Below phreatic surface ---')
    print(f'Max kr: {np.max(kr[below_mask]):.3e}')
    print(f'Mean kr: {np.mean(kr[below_mask]):.3e}')
    print(f'Min kr: {np.min(kr[below_mask]):.3e}')

    # Plot
    plt.figure(figsize=(10, 8))
    plt.scatter(coords[below_mask, 0], coords[below_mask, 1], c=kr[below_mask], cmap='Blues', label='Below Phreatic', s=60)
    plt.scatter(coords[above_mask, 0], coords[above_mask, 1], c=kr[above_mask], cmap='Reds', label='Above Phreatic', s=60, marker='^')
    plt.colorbar(label='kr')
    plt.legend()
    plt.title(title)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid(True)
    plt.axis('equal')
    plt.show()

def debug_kr_elements(coords, elements, head, kr_elem, title='Element kr vs Phreatic Surface'):
    """
    Plot and summarize kr values at the element level, classified by phreatic surface.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    # Compute average head and average elevation for each element
    avg_head = np.mean(head[elements], axis=1)
    avg_elev = np.mean(coords[elements, 1], axis=1)
    above_mask = avg_head < avg_elev
    below_mask = ~above_mask

    print(f'Elements above phreatic surface: {np.sum(above_mask)}')
    print(f'Elements below phreatic surface: {np.sum(below_mask)}')
    print('--- Above phreatic surface ---')
    print(f'Max kr: {np.max(kr_elem[above_mask]):.3e}')
    print(f'Mean kr: {np.mean(kr_elem[above_mask]):.3e}')
    print(f'Min kr: {np.min(kr_elem[above_mask]):.3e}')
    print('--- Below phreatic surface ---')
    print(f'Max kr: {np.max(kr_elem[below_mask]):.3e}')
    print(f'Mean kr: {np.mean(kr_elem[below_mask]):.3e}')
    print(f'Min kr: {np.min(kr_elem[below_mask]):.3e}')

    # Plot
    verts = [coords[tri] for tri in elements]
    fig, ax = plt.subplots(figsize=(10, 8))
    pc = PolyCollection(verts, array=kr_elem, cmap='viridis', edgecolor='k')
    ax.add_collection(pc)
    ax.autoscale()
    plt.colorbar(pc, ax=ax, label='kr')
    plt.title(title)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid(True)
    plt.show()

def debug_ke_elements(coords, elements, head, ke_list, title='Element ke vs Phreatic Surface'):
    """
    Print and compare the norm of ke for elements above vs below the phreatic surface.
    ke_list: list or array of 3x3 element stiffness matrices
    """
    import numpy as np
    # Compute average head and average elevation for each element
    avg_head = np.mean(head[elements], axis=1)
    avg_elev = np.mean(coords[elements, 1], axis=1)
    above_mask = avg_head < avg_elev
    below_mask = ~above_mask

    # Compute Frobenius norm of each ke
    ke_norms = np.array([np.linalg.norm(ke) if ke is not None else np.nan for ke in ke_list])

    print(f'Elements above phreatic surface: {np.sum(above_mask)}')
    print(f'Elements below phreatic surface: {np.sum(below_mask)}')
    print('--- Above phreatic surface ---')
    print(f'Max ke norm: {np.nanmax(ke_norms[above_mask]):.3e}')
    print(f'Mean ke norm: {np.nanmean(ke_norms[above_mask]):.3e}')
    print(f'Min ke norm: {np.nanmin(ke_norms[above_mask]):.3e}')
    print('--- Below phreatic surface ---')
    print(f'Max ke norm: {np.nanmax(ke_norms[below_mask]):.3e}')
    print(f'Mean ke norm: {np.nanmean(ke_norms[below_mask]):.3e}')
    print(f'Min ke norm: {np.nanmin(ke_norms[below_mask]):.3e}')