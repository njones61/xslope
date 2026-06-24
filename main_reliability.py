from xslope.fileio import load_slope_data
from xslope.plot import plot_reliability_results
from xslope.advanced import reliability

# Load slope data
slope_data = load_slope_data("docs/inputs/slope/input_template_reliability6.xlsx")

# Run reliability analysis. `method` can be any LEM method: 'oms', 'bishop',
# 'janbu', 'corps', 'lowe', 'spencer', or 'mprice'.
print("Running reliability analysis with Spencer method...")
success, result = reliability(slope_data, method='spencer', rapid=False, circular=True, debug_level=1)

if success:
    print("\nReliability analysis completed successfully!")
    
    # Plot the results
    plot_reliability_results(slope_data, result)
    
else:
    print(f"Reliability analysis failed: {result}")