from slice import generate_slices
from fileio import load_globals

# Load test data using the same method as main.py
data = load_globals("docs/input_template_dam_copy.xlsx")

# Print profile_lines and materials info for debugging
print('Num profile_lines:', len(data['profile_lines']))
print('Num materials:', len(data['materials']))
for i, mat in enumerate(data['materials']):
    print(f'Material {i+1}:', mat)

# Print profile_lines coordinates
print('\nProfile lines:')
for i, profile in enumerate(data['profile_lines']):
    print(f'Profile {i}: {profile}')

# Get the circle data
circle = data['circles'][0] if data['circular'] else None

# Run generate_slices
success, result = generate_slices(data, circle=circle, debug=True)

print('Success:', success)
if success:
    df, surface = result
    # Find all h columns
    h_cols = [col for col in df.columns if col.startswith('h') and col[1:].isdigit()]
    h_cols = sorted(h_cols, key=lambda x: int(x[1:]))
    print('\nFirst few slices:')
    print(df[['slice #', 'x_c'] + h_cols].head())
    
    print('\nAll slices h values:')
    for i, row in df.iterrows():
        h_vals = ', '.join([f"{col}={row[col]:.3f}" for col in h_cols])
        print(f"Slice {row['slice #']}: {h_vals}")
else:
    print('Error:', result) 