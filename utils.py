from shapely.geometry import LineString, Point

def build_ground_surface(profile_lines):
    # Step 1: Gather all points and sort by x, then by descending y
    all_points = sorted(set(pt for line in profile_lines for pt in line), key=lambda p: (p[0], -p[1]))

    # Step 2: Keep only the highest y for each x
    top_candidates = {}
    for x, y in all_points:
        if x not in top_candidates or y > top_candidates[x]:
            top_candidates[x] = y

    # Step 3: For each top candidate, check that it's higher than all other lines
    top_surface_points = []
    for x, y in sorted(top_candidates.items()):
        keep = True
        for other_line in profile_lines:
            line = LineString(other_line)
            if line.length == 0:
                continue
            proj = line.project(Point(x, 0))
            if proj == 0 or proj == line.length:
                continue  # avoid extrapolation
            ipt = line.interpolate(proj)
            if ipt.y > y + 1e-6:
                keep = False
                break
        if keep:
            top_surface_points.append((x, y))

    return LineString(top_surface_points) if len(top_surface_points) >= 2 else LineString([])
