from fileio import load_slope_data
from plot import plot_reliability_results
from advanced import reliability
from solve import spencer, bishop

# Load slope data
slope_data = load_slope_data("inputs/slope/input_template_reliability4.xlsx")

# Run reliability analysis with Spencer method
print("Running reliability analysis with Spencer method...")
success, result = reliability(slope_data, method=spencer, rapid=False, circular=True, debug_level=1)

if success:
    print("\nReliability analysis completed successfully!")
    
    # Plot the results
    plot_reliability_results(slope_data, result)
    
else:
    print(f"Reliability analysis failed: {result}")