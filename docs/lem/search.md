# Automated Search Algorithms

The xslope package provides sophisticated automated search algorithms for locating critical failure surfaces in slope stability analysis. These algorithms systematically explore the design space to find the surface geometry that yields the minimum factor of safety, which represents the most dangerous potential failure mode. The search functionality supports both circular and non-circular failure surfaces, each employing distinct optimization strategies tailored to the geometric constraints of the failure mechanism.

## Circular Failure Surfaces

The circular search algorithm (`circular_search()` in xslope/search.py) implements a global optimization approach using adaptive grid refinement to locate the critical circular failure surface. The algorithm begins with one or more user-defined starting circles that serve as initial guesses for the optimization process. Each starting circle is characterized by its center coordinates (Xo, Yo), depth, and radius. The initial guess may be approximate; the algorithm refines it toward the critical surface.

The circular search uses a two-stage optimization that operates on both the location of the circle center and the depth of the failure surface. When evaluating a particular center location (x, y), the algorithm first performs a nested depth optimization: for any given center point, there is a depth that minimizes the factor of safety. The depth optimizer uses a three-point bracketing approach, evaluating the factor of safety at the current best depth and at depths one step above and below it. By iteratively selecting the depth that produces the lowest factor of safety and shrinking the step size by a factor of 0.25, the algorithm converges to the locally optimal depth for that center location. This depth optimization continues until the step size falls below a tolerance threshold, typically set to 1% of the current step size.

Once the algorithm can efficiently evaluate any point in the (x, y) plane by optimizing its depth, it employs a nine-point grid search strategy to locate the optimal center coordinates. Starting from the best point found among all initial starting circles, the algorithm constructs a 3×3 grid centered on the current best location. The grid spacing is initially set to 15% of the circle radius, providing a balance between exploration and computational efficiency. The algorithm evaluates the factor of safety at all nine points on this grid, utilizing a caching mechanism to avoid redundant calculations when grid points overlap between iterations. 

![alt text](images/autosearch_grid1.png)

After identifying which of the nine points yields the lowest factor of safety, the algorithm updates its center location if an improvement was found. For example, for the grid shown above, suppose the upper right point had the lowest FS of the 9 grid points. The next grid would then shift so that this is the center point:

![alt text](images/autosearch_grid2.png)

The new points on this grid are then evaluated with the circles at various depths to get a new FS value for each point. Then if the grid point on the middle right had the lowest FS, the grid would shift again as follows:

![alt text](images/autosearch_grid3.png)

Again, the new points on this grid are evaluated for minimum FS. If the best point on the grid is the center point itself (meaning no improvement was discovered), the algorithm shrinks the grid size by multiplying it by the shrink factor (default 0.5), effectively zooming in on the current location to search at a finer resolution as follows:

![alt text](images/autosearch_grid4.png)

Each of the new points on this grid are then evaluated. This iterative process of grid evaluation and refinement continues until the **factor of safety stops changing**, rather than until the grid reaches some fixed spacing. Specifically, the search converges once two successive refinement levels each lower the best factor of safety by less than `fs_tol` (default `5e-4`, so the reported FS is stable to three decimal places), provided the center grid has already refined below a small fraction of the slope height (`min_grid_frac`, default 1%). That second condition prevents a coarse grid from stopping early on a plateau before it has resolved the true minimum. A geometric grid floor (`tol`) and an iteration cap (`max_iter`) remain only as backstops. Keying convergence on the factor of safety rather than an absolute grid length makes the stopping point **scale-invariant** — a fixed length threshold such as "grid < 0.01" would mean something entirely different on a 20 ft slope than on a 500 ft slope. The algorithm tracks its progression through the search space, recording each jump to a new best location along with the corresponding factor of safety. Diagnostic output reports the center coordinates, factor of safety, and grid size at each iteration. The search typically converges within 10-20 iterations, depending on the complexity of the slope geometry and the quality of the initial starting guess.

Throughout the search process, the algorithm maintains a comprehensive cache (`fs_cache`) that stores every evaluated circle configuration along with its computed factor of safety and associated solution details. This cache serves multiple purposes: it prevents redundant calculations, enables post-analysis visualization of the entire search space, and provides a ranked list of all evaluated failure surfaces. The cache is implemented as a dictionary keyed by (x, y) coordinates, with each entry storing the circle center coordinates, depth, factor of safety, slice dataframe, failure surface geometry, and the detailed solver results.

At the end of the process, the search results can be displayed using the **plot_circular_search_results** function as follows:

![alt text](../seep/images/seep_slope_lem_search.png)

Note that all circles tested in the search are displayed along with all grid points evaluated using the 9-point search algorithm. Furthermore, the sequence of center points for the moving 9-point grids are connected with a series of arrows showing the "search path".

### Searches with Several Starting Circles {#several-starting-circles-are-not-several-searches}

The circles sheet may hold any number of starting circles, but a seeded search
(`seed='circles'`, the default) does not refine them all. Each starting circle is
screened once — a single nine-point grid at its own launch spacing — the circles are
ranked by that one coarse score, and the adaptive refinement then runs from the
best-scoring circle alone. The other circles' families are never explored, and
nothing in the output says so: the run reports convergence, one factor of safety,
and one critical surface. Two circles entered to check two mechanisms, such as an
upstream and a downstream surface on a dam, therefore report one mechanism, and
the one that is kept is chosen at the coarse resolution of the screening grid: a
family that screens a few thousandths higher but would refine much lower is dropped
before refinement begins. Enter a single starting circle, on the face known to
control, when the run is meant to examine one specific mechanism;
turn on grid seeding, described next, when it is meant to find the critical surface
anywhere in the model, since that mode refines up to four competing families with
the sheet's own circles among them.

A starting circle that cannot be sliced — most often one whose crossing with the
ground lies above the circle's own center, so the arc would have to rise above
its center on that side — is not searched from. The log reports this
(`starting circle 1 … cannot be built … searching from the launch grid`) and the
search continues from the grid of centers around the circle; the model
checks report the same circle as `surface.circle_daylights_above_center` before
the run.

### Grid Seeding (Global Search)

The adaptive search described above is a **local** optimizer: it refines whatever neighborhood its starting circles put it in, and it will do so even when a far better minimum exists elsewhere. On an embankment over a layered foundation, a shallow circle in the fill and a deep circle along the base of the foundation are two *competing families* of failure surface, separated by a ridge of higher factors of safety that the local search will not cross. Started from a single circle in the wrong family, the search converges cleanly — reporting convergence, a stable factor of safety, and a plausible-looking surface — at a value that can be 20% or more too high. The output gives no warning of this.

Grid seeding (`seed='grid'`) addresses this by adding a global stage in front of the local search:

```python
fs_cache, converged, search_path, circle_cache = circular_search(
    slope_data, 'spencer', seed='grid')
```

The global stage is the classical **grid-and-tangent** sweep (the model used by Slide2 and SLOPE/W). Circle centers are laid out on a coarse grid over a box derived from the slope geometry — spanning the slope horizontally and sitting above the crest at heights scaled to the full crest-to-base depth of the model, which brackets where the center of any circle that cuts the slope can be. For each center, one trial circle is made *tangent* to each of a series of elevations between the bottom of the model and mid-slope. Sweeping tangent elevations rather than depths separates competing families: the shallow face circle and the deep base-tangent circle occur at similar centers but different tangent lines, so each family appears as its own row of the sweep. Every trial circle (typically ~300) is solved at a reduced slice count, with the same degenerate-surface guards as the main search.

The best circle from each surviving tangent family — every family within 25% of the best factor of safety found, up to four — then seeds the adaptive search, together with any circles from the input file, and **each family is refined independently**. Families can differ by only a few percent at the coarse stage yet refine to very different minima, so refining only the best coarse result would bring back the same problem at the next stage. The shared evaluation cache makes the extra refinements cheap; the whole grid-seeded search typically costs a few times a plain seeded search (seconds to tens of seconds).

Because the seeds come from the geometry, grid seeding also works with an **empty circles sheet** — the search needs no starting guess at all.

The global search differs from the seeded search in two ways:

- **It reports the most critical surface anywhere in the model.** If the model contains a genuinely marginal minor feature — say a steep fill face standing near its limit — the global minimum is there, not on the main slope, and that is what the search returns. This is correct behavior, but it differs from a benchmark or a design check that targets one specific mechanism. To examine a particular mechanism, use the default seeded mode with circles placed in that mechanism's family.
- **Surficial skin slides are filtered — in grid mode only.** On a steep face with little cohesion, a very thin surface along the face can give a lower factor of safety than every real mechanism (a raveling mode, not a slope failure), so grid mode rejects trial surfaces whose maximum thickness is under 5% of the slope height. The filter deliberately does **not** apply to the default seeded search, because a thin surface can also be the correct answer — a submerged cohesionless slope fails as an infinite-slope skin at essentially zero depth, and a seeded search must be able to follow your circles there. If the critical surface of a problem is that kind of skin, use seeded mode.

On the James Bay dyke verification problem ([VP75](../verification/rocscience.md#vp75)), a single plausible starting circle leaves the local search at FS = 1.745; grid seeding finds FS = 1.421 from the same file with the circles sheet ignored entirely — marginally better than the 1.425 found from three hand-picked seeds, and closest to Slide2's published 1.415.

In XSLOPE Studio, grid seeding is the **Grid search** toggle in the Run LEM dialog, available when the analysis is an automated circular search.

### Search-Window Limits

By default, `circular_search()` places no limits on where a trial circle may lie. Four optional keyword arguments — `center_box`, `entry_range`, `exit_range`, and `tangent_depth` — narrow that exploration to a chosen region instead, the LEM analog of Slide2's slip-center / entry-and-exit search limits and RS2's SSR Polygon Search Area. All four default to `None` (unconstrained — the ordinary behavior of `circular_search()`); passing any of them confines the search to a region so it settles on a chosen **local** minimum rather than the global one.

- `center_box=(x1, y1, x2, y2)` restricts candidate circle **centers** to a rectangle (corners in any order). Grid points outside the box are dropped outright, so refinement can never walk out of it; a starting circle whose center lies outside is clamped in only for its own launch point, not reported as a result.
- `entry_range=(xa, xb)` / `exit_range=(xc, xd)` restrict where the failure surface crosses the ground surface, in x. Entry is the crest-side (higher-ground) endpoint and exit the toe-side one, independent of slope facing. A trial whose endpoints fall outside its range is **rejected** (scored at `fs_fail`), never clamped, so the reported minimum lies within the window.
- `tangent_depth=(ymin, ymax)` restricts the elevation of the circle's lowest (tangent) point to a band; out-of-band trials are rejected the same way.

All four constraints apply in both the adaptive 9-point refinement and the `seed='grid'` sweep, and a malformed bound raises `ValueError`.

```python
fs_cache, converged, search_path, circle_cache = circular_search(
    slope_data, 'spencer',
    entry_range=(42, 54),
    exit_range=(23, 32),
    tangent_depth=(16, 22),
)
```

Use these limits when a slope has several competing local minima and the question at hand is about a specific one of them, not whichever is lowest. [RS2-61](../verification/rs2.md#rs2-61) is a worked example: the same homogeneous-slope geometry has a global minimum (Spencer 1.338, matching Slide2's 1.336) and a separate upper-face local minimum (Spencer 1.437, matching Slide2's 1.443) that an unconstrained search never reports, because it always finds the lower global value first. An `entry_range` / `exit_range` / `tangent_depth` window isolates the upper-face family instead.

!!! note "API only"
    These search-window limits are available at the `circular_search()` / test-tag level only — there is no dialog control for them in XSLOPE Studio. Set them from the API, or via a `circular_search` test tag (`center_box`, `entry_range`, `exit_range`, `tangent_depth`, each a `;`-separated list of numbers).

## Non-Circular Failure Surfaces

The non-circular search algorithm (`noncircular_search()` in xslope/search.py) takes a fundamentally different approach from the circular search, as it must optimize the positions of multiple control points that define an arbitrary failure surface rather than just optimizing a center and depth. The algorithm begins with a user-defined initial failure surface specified as a sequence of points in the slope_data['non_circ'] array. Each point in this sequence is characterized by its x and y coordinates along with a movement constraint that controls how the point is allowed to move during the optimization process.

The movement constraints are critical to ensuring that the optimized failure surface remains physically realistic and geometrically valid. Points can be designated as "Fixed" (remaining stationary throughout the search), "Horiz" (allowed to move only horizontally), or "Free" (able to move in any direction that reduces the factor of safety). The endpoints of the failure surface have special handling: a `Fixed` end vertex is pinned like any other fixed point, and a movable one is constrained to the ground surface itself — the search steps it by arc length *along* the ground polyline, so an end walks the surface it daylights on, including up or down a vertical face such as a wall. `Free` and `Horiz` are equivalent on an end vertex, since along-the-ground is the only admissible direction there. This ensures that the failure surface always properly intersects the slope at both ends.

The optimization algorithm employs a coordinate descent strategy, systematically attempting to move each point to reduce the overall factor of safety. During each iteration, the algorithm cycles through all points in the failure surface, and for each movable point, it tests movements in both the positive and negative directions. For points constrained to horizontal movement, these directions are simply left and right displacements. For free points, however, the algorithm implements a more sophisticated movement strategy: it moves each free point perpendicular to the local tangent of the failure surface. This perpendicular movement is calculated by finding the tangent vector between the point's neighbors, rotating it 90 degrees, and normalizing it to the current movement distance. This approach naturally tends to smooth out irregularities in the surface while allowing it to deform toward lower factors of safety.

The movement distance parameter controls the step size of these trial movements, starting at a user-specified initial value (typically 4.0 length units). After attempting to move each point in both directions, the algorithm checks whether any improvement in the factor of safety was achieved. If the entire pass through all points fails to find any improvement, or if the change in factor of safety falls below the tolerance threshold (default 0.001), the movement distance is reduced by multiplying it by the shrink factor (default 0.8). This gradual reduction in step size allows the algorithm to make large exploratory moves early in the search while refining the solution with increasingly precise adjustments as it approaches convergence.

Geometric validity is strictly enforced through a series of constraints checked by the `move_point()` function. The x-coordinates of points must maintain their sequential ordering—each point must remain to the right of its left neighbor and to the left of its right neighbor. Middle points (non-endpoints) must remain below the ground surface, which is verified by interpolating the ground surface elevation at the point's x-coordinate and ensuring the y-coordinate doesn't exceed it. Additionally, all points must remain above the maximum depth constraint specified in the slope geometry. If a proposed movement would violate any of these constraints, it is rejected and the point remains in its current position.

A further **base-angle limit** rejects any trial surface whose base inclination exceeds `max_base_angle` (default 65°, the active-wedge angle 45° + φ/2 for φ ≈ 40°, near the steep limit of realistic soils). Without this cap the coordinate-descent search tends to drive the points into a near-vertical base segment running up to the toe — a geometrically unrealistic surface that the **force-equilibrium** methods, Lowe-Karafiath above all, score as a spurious low minimum, because their interslice inclination is tied to the base slope. The complete-equilibrium methods are not affected: on the weak-layer model of [Tutorial LEM-5](../tutorials/lem05_weak_layer_noncircular.md) Spencer and Morgenstern-Price return the same factor of safety whether the cap stands at 65° or is opened to 89°, and the surface each settles on is nowhere near either limit, while Lowe-Karafiath's minimum falls step by step and its surface sits on whatever limit is set. Capping the base angle keeps the search on physically realistic surfaces; trial surfaces that exceed the cap are treated as inadmissible (assigned an infinite factor of safety) so the search never settles on them.

Convergence is determined using an AND logic that requires both the change in factor of safety between iterations to fall below the specified tolerance (default 0.001) AND the movement distance to shrink below the movement tolerance (default 0.1). This dual criterion ensures that the algorithm has not only found a local minimum in terms of factor of safety but has also refined the surface position to a satisfactory precision. The maximum iteration limit (default 100) provides a safeguard against excessive computation time, though most searches converge well before reaching this limit.

Once again, XSLOPE preserves the search results and all of the non-circular results tested can be plotted with the **plot_noncircular_search_results** function.

![alt text](sample_images/noncircular_search_results.png)

## Search Results Structure

The circular search returns a four-element tuple `(fs_cache, converged, search_path, circle_cache)`, and the non-circular search returns a three-element tuple `(fs_cache, converged, search_path)`. These return structures provide comprehensive information about the optimization process, supporting both immediate access to the critical failure surface and detailed post-processing analysis of the entire search trajectory.

The `fs_cache` is returned as a sorted list of dictionaries, ordered by factor of safety from lowest to highest. Each dictionary in this list represents one evaluated failure surface configuration and contains several key fields. The 'FS' field stores the computed factor of safety, which serves as the primary metric for ranking surfaces. The 'slices' field contains the complete pandas DataFrame of slice properties that was used in the limit equilibrium calculation, including slice widths, base inclinations, weights, normal forces, and shear strengths. The 'failure_surface' field holds the shapely LineString geometry representing the exact path of the failure surface through the slope. The 'solver_result' field stores the detailed output from the limit equilibrium solver method (such as Bishop, Spencer, or Janbu), which may include additional information like interslice force angles, iteration counts, and convergence metrics.

For circular searches, each cache entry additionally includes the 'Xo', 'Yo', and 'Depth' fields that fully specify the circle geometry. For non-circular searches, the 'points' field contains the array of control point coordinates. This structure means that the critical failure surface—the one with the minimum factor of safety—is always the first element in the returned `fs_cache[0]`, making it trivial to extract the most important result.

The `converged` boolean flag indicates whether the search algorithm successfully converged according to its termination criteria or whether it hit the maximum iteration limit. A value of `True` means the algorithm found a local minimum and refined it to within the specified tolerances, while `False` suggests that either more iterations are needed or the algorithm encountered numerical difficulties. Users should check this flag and potentially adjust search parameters if convergence was not achieved.

The `search_path` list documents the progression of the search algorithm through the design space, recording each significant step where a better solution was found. For circular searches, each element in the search path is a dictionary with 'x', 'y', and 'FS' fields showing where the circle center moved and what factor of safety was achieved. For non-circular searches, the path tracking is more complex, as it needs to capture the evolution of the entire multi-point surface. This search path data enables visualization of the optimization trajectory and helps diagnose potential issues like premature convergence to local minima or oscillatory behavior.

The circular search additionally returns a `circle_cache` containing every circle evaluated during the search (not just the improving steps), which is used to plot all tested circles alongside the critical surface.

### Trial Surfaces the Method Could Not Solve

A trial surface can be geometrically admissible and still have no solution by the selected method. Every method assumes something about the interslice forces, and each assumption has surfaces it cannot satisfy. Such a trial is scored the same as one rejected on geometry and dropped from the ranking, so the reported minimum is the minimum of the surfaces the method could solve.

The circular search counts these trials. When the method could not solve one or more admissible trials, the search prints one line giving the count and the reason. This line is from Spencer on the embankment of [Tutorial LEM-1](../tutorials/lem01_simple_embankment.md), where the material is undrained (φ = 0):

```
[⚠️ unsolved trials] Spencer could not solve 56 of 211 trial surfaces (on 56, no
interslice force inclination satisfies both force and moment equilibrium)
```

A search in which the method solved every admissible trial prints nothing extra. Pass a dictionary as `unsolved_out=` to receive the same counts (`attempted`, `unsolved`, `no_admissible_solution`, `not_converged`, `inadmissible`, `unclassified`, and the printed `sentence`); `run_lem_analysis` carries them on `search["unsolved"]`, which is `None` for a non-circular search, whose scoring does not separate a solver failure from a geometric rejection.

The reasons in parentheses are of three kinds.

- **No solution.** The method's equations have no acceptable solution on the surface. Spencer prints "no interslice force inclination satisfies both force and moment equilibrium": no single inclination within the allowable range closes both equations. Roots can exist outside that range, but there they reverse the normal force on the base of some slices, so they are not considered. Morgenstern-Price prints "no value of λ from -1.5 to 1.5 satisfies both force and moment equilibrium": its force and moment factors of safety do not cross anywhere in that range of the interslice force scale factor λ. The Corps of Engineers and Lowe & Karafiath methods print "no factor of safety gives a physically reasonable solution to force equilibrium": either no factor of safety from 0.05 to 10 closes force equilibrium, or every one that does fails the tests listed on the [force-equilibrium method page](force_eq.md).
- **Failed to converge.** The iteration did not reach a solution that may well exist.
- **Solution rejected as physically unreasonable.** The method found a solution and rejected it. For Spencer's method a solution is rejected if any of three tests fails:
    1. More than half of the slice bases are in tension (negative effective normal force).
    2. The interslice force inclination lies outside its allowable range. The range is the one in which every slice keeps a positive base normal force; it corresponds to an inclination between the slope of the ground surface and the slope of the slip surface.
    3. An interslice force is larger than the total load on the sliding mass: its weight plus the distributed loads, line loads, seismic force and tension-crack water force applied to it.

    A solution with a factor of safety of zero or less is also rejected. Morgenstern-Price applies the first test only. The Ordinary Method of Slices, Bishop's method and the Janbu method apply the first test, and Bishop's and Janbu's methods also reject a solution in which m<sub>α</sub> falls below 0.05 on any slice.

The breakdown appears only when XSLOPE has classified every failure message the method returned. Otherwise the line gives the count alone. Interslice tension does not cause a trial to be rejected; it is [reported instead](overview.md#interslice-tension).

Two changes can help when many trials go unsolved. Another method's assumption may be satisfiable where this one's is not, so running a second method shows how much of the model the first method's answer covers. A surface with no solution because the soil near the crest is in tension becomes solvable once the model includes a tension crack, since the crack removes the slices whose cohesion the assumption cannot balance.

### Tension at the Crest

When the reported solution has substantial tension on the base of the slices at the crest end of the surface, XSLOPE adds a note to the solution's warnings. If the model has no tension crack, the note suggests one:

```
Tension on the base of 3 slices near the crest. Consider adding a tension crack.
```

The crest end is the end of the slip surface where it meets the higher ground. Counting from that end, the run is the slices that are in base tension (negative effective normal force) and have cohesion; it stops at the first slice that is not. The note appears only when all of the following hold:

- The run has at least 2 slices.
- The worst base tension in the run, as a stress (effective normal force divided by base length), is at least half the cohesion of that slice.
- No reinforcement line or pile acts on a slice of the run or on the slice next to it.

Slices whose base is shorter than 0.5% of the length of the slip surface are passed over in the count. If the model already has a tension crack and the run meets these conditions, the note gives the first sentence only. Smaller tension at the crest, tension elsewhere on the surface, and tension on cohesionless slices draw no note of this kind; the method's own notes report base tension on cohesionless slices. The note appears with the other warnings: in the run output, where the warnings print before the factor of safety, in the amber strip in Studio, in the report and in the assistant's account of the run. It does not change the factor of safety or which surface is reported. [Tutorial LEM-1](../tutorials/lem01_simple_embankment.md#the-fix-is-in-the-ground-not-the-settings) shows the remedy.

## Visualization of Search Results

The xslope package provides specialized plotting functions that transform the numerical search results into intuitive visual representations of the failure surface exploration process. These functions use matplotlib to create publication-quality figures that overlay the discovered failure surfaces onto the slope geometry.

The `plot_circular_search_results()` function (in xslope/plot.py) creates a comprehensive visualization that begins by plotting all the fundamental slope features—the ground surface profile, subsurface stratigraphy boundaries, maximum depth line, piezometric surface, distributed loads, and tension crack geometry. This establishes the physical context within which the failure surfaces were searched. A piezometric line is drawn only where the analysis reads it — as a material's pore-pressure source (`u = piezo`), as the water table a `gamma_sat` weight split is measured from, or as the sheet that states the pool loading the slope — so a line the analysis does not use does not appear on the figure. The function then overlays every circular failure surface stored in the `fs_cache`, plotting them as curved lines. The critical surface (with the minimum factor of safety) is rendered in a bold red line with width 2, while all other explored surfaces are shown as semi-transparent gray lines with reduced width. This visual hierarchy immediately draws attention to the most dangerous surface while still revealing the breadth of the search space that was explored.

The circle centers are marked with small points, and if the `search_path` was provided, the function plots arrows showing how the optimization algorithm moved from one center location to another. These arrows create a visual trace of the search trajectory, making it easy to see whether the algorithm made large jumps early on before converging with small refinements, or whether it struggled with many small steps. The plot title displays the critical factor of safety value, providing immediate quantitative feedback. The resulting figure uses equal aspect ratio to prevent geometric distortion and includes a legend to identify the different surface types.

Similarly, `plot_noncircular_search_results()` (in xslope/plot.py) follows the same pattern for non-circular surfaces, plotting all cached failure surfaces with the critical one highlighted. For non-circular searches, the search path visualization is more elaborate—instead of simple arrows between center points, it can show arrows for each control point that moved between iterations, revealing how different parts of the failure surface evolved during optimization. This detailed movement tracking helps users understand which portions of the surface were most active in the search and whether the final surface has smooth, natural-looking geometry or potentially unrealistic kinks that might indicate convergence issues.

Both plotting functions accept optional parameters for figure size, whether to save the plot as a PNG file, and the DPI resolution for saved figures. The `highlight_fs` parameter can be set to `False` to give all surfaces equal visual weight, which is useful when examining the entire population of explored surfaces rather than focusing solely on the critical case.

## Example Usage

The following example demonstrates a complete workflow for performing an automated circular failure surface search and visualizing the results:

```python
from xslope.fileio import load_slope_data
from xslope.search import circular_search
from xslope.plot import plot_circular_search_results

# Load slope geometry and material properties from Excel template
slope_data = load_slope_data("inputs/slope/input_template_reliability6.xlsx")

# Perform circular search using Spencer's method
# The 'circles' in slope_data provide the starting guess(es)
fs_cache, converged, search_path, circle_cache = circular_search(
    slope_data,
    method_name='spencer',
    diagnostic=False,    # Set True to see detailed iteration output
    num_slices=40,       # Slices per trial surface (default 40)
    fs_tol=5e-4,         # FS convergence tolerance (stable to ~3 decimals)
    min_grid_frac=0.01,  # Min center-grid resolution (fraction of slope height)
    max_iter=50,         # Maximum iterations (backstop)
    shrink_factor=0.5,   # Grid reduction factor
    tol=0.01,            # Geometric grid floor (backstop only)
)

# Check if search converged
if converged:
    print(f"Search converged successfully")
else:
    print(f"Warning: Search did not converge within max iterations")

# Extract critical failure surface
critical_surface = fs_cache[0]
critical_fs = critical_surface['FS']
critical_center = (critical_surface['Xo'], critical_surface['Yo'])
critical_depth = critical_surface['Depth']

print(f"Critical FS = {critical_fs:.4f}")
print(f"Circle center at ({critical_center[0]:.2f}, {critical_center[1]:.2f})")
print(f"Failure depth = {critical_depth:.2f}")

# Visualize all explored surfaces with critical one highlighted
plot_circular_search_results(slope_data, fs_cache, search_path, save_png=True)

# Access detailed solution information
slices_df = critical_surface['slices']
solver_result = critical_surface['solver_result']
print(f"\nSolver converged in {solver_result.get('iterations', 'N/A')} iterations")
```

For non-circular failure surface searches, the workflow is nearly identical but uses `noncircular_search()` and `plot_noncircular_search_results()`:

```python
from xslope.fileio import load_slope_data
from xslope.search import noncircular_search
from xslope.plot import plot_noncircular_search_results

# Load slope data (must include 'non_circ' specification)
slope_data = load_slope_data("inputs/slope/input_template_noncircular.xlsx")

# Perform non-circular search using Lowe-Karafiath method
fs_cache, converged, search_path = noncircular_search(
    slope_data,
    method_name='lowe',
    diagnostic=True,           # Print iteration details
    movement_distance=4.0,     # Initial step size
    shrink_factor=0.8,         # Step reduction factor
    fs_tol=0.001,             # FS convergence tolerance
    max_iter=100,             # Maximum iterations
    move_tol=0.1              # Movement distance tolerance
)

# Extract and display critical surface
critical_surface = fs_cache[0]
critical_points = critical_surface['points']
critical_fs = critical_surface['FS']

print(f"Critical FS = {critical_fs:.4f}")
print(f"Failure surface defined by {len(critical_points)} points:")
for i, (x, y) in enumerate(critical_points):
    print(f"  Point {i}: ({x:.2f}, {y:.2f})")

# Visualize the search results
plot_noncircular_search_results(slope_data, fs_cache, search_path)
```

Both search algorithms can also be combined with rapid drawdown analysis by setting the `rapid=True` parameter, which modifies the solution method to account for transient pore pressure conditions during reservoir drawdown scenarios. The search algorithms pass the appropriate parameters to the underlying limit equilibrium solvers automatically, so the search itself needs no other change.