import numpy as np
import pandas as pd
from fileio import load_slope_data
from slice import generate_slices

def numerical_derivative(func, x, h=1e-6):
    """Compute numerical derivative using central difference."""
    return (func(x + h) - func(x - h)) / (2 * h)

def numerical_second_derivative(func, x, h=1e-6):
    """Compute numerical second derivative."""
    return (func(x + h) - 2 * func(x) + func(x - h)) / (h**2)

# Load test data
slope_data = load_slope_data('inputs/slope/input_template_lface4.xlsx')
success, slice_df = generate_slices(slope_data)
if not success:
    print("Failed to generate slices")
    exit(1)

# Extract Spencer method variables (simplified)
alpha = np.radians(slice_df['alpha'].values)
phi = np.radians(slice_df['phi'].values)
c = slice_df['c'].values
dx = slice_df['dx'].values
dl = slice_df['dl'].values
W = slice_df['w'].values
u = slice_df['u'].values
x_c = slice_df['x_c'].values
y_cb = slice_df['y_cb'].values
P = slice_df['dload'].values
beta = np.radians(slice_df['beta'].values)
kw = slice_df['kw'].values
V = slice_df['t'].values
y_v = slice_df['y_t'].values
R = slice_df['p'].values

# For now, we assume that reinforcement is flexible and therefore is parallel to the failure surface
psi = alpha
y_r = y_cb
x_r = x_c

x_p = slice_df['d_x'].values
y_p = slice_df['d_y'].values
y_k = slice_df['y_cg'].values
x_b = x_c
y_b = y_cb

tan_p = np.tan(phi)

# Check if right facing
y_ct = slice_df['y_ct'].values
right_facing = (y_ct[0] > y_ct[-1])

if right_facing:
    alpha = -alpha
    beta = -beta
    psi = -psi
    R = -R
    c = -c
    kw = -kw
    tan_p = -tan_p

# Pre-compute trigonometric functions
cos_a = np.cos(alpha)
sin_a = np.sin(alpha)
cos_b = np.cos(beta)
sin_b = np.sin(beta)
sin_psi = np.sin(psi)
cos_psi = np.cos(psi)

Fh = -kw - V + P * sin_b + R * cos_psi
Fv = -W - P * cos_b + R * sin_psi
Mo = -P * sin_b * (y_p - y_b) - P * cos_b * (x_p - x_b) + kw * (y_k - y_b) + V * (y_v - y_b) - R * cos_psi * (y_r - y_b) + R * sin_psi * (x_r - x_b)

# Test point
F_test = 1.5
theta_test = np.radians(8.0)

def compute_Q(F, theta_rad):
    """Compute Q for given F and theta."""
    ma = 1 / (np.cos(alpha - theta_rad) + np.sin(alpha - theta_rad) * tan_p / F)
    Q = (-Fv * sin_a - Fh * cos_a - (c / F) * dl + (Fv * cos_a - Fh * sin_a + u * dl) * tan_p / F) * ma
    return Q

def compute_y_q(F, theta_rad):
    """Compute y_q for given F and theta."""
    Q = compute_Q(F, theta_rad)
    y_q = y_b + Mo / (Q * np.cos(theta_rad))
    return y_q

def compute_R1(F, theta_rad):
    """Compute R1 residual."""
    Q = compute_Q(F, theta_rad)
    return np.sum(Q)

def compute_R2(F, theta_rad):
    """Compute R2 residual."""
    Q = compute_Q(F, theta_rad)
    y_q = compute_y_q(F, theta_rad)
    return np.sum(Q * (x_b * np.sin(theta_rad) - y_q * np.cos(theta_rad)))

# Test numerical vs analytical derivatives
print("Testing derivatives at F=1.5, theta=8°")
print()

# Test dR1/dF
dR1_dF_numerical = numerical_derivative(lambda F: compute_R1(F, theta_test), F_test)
print(f"dR1/dF numerical: {dR1_dF_numerical:.6e}")

# Test dR1/dtheta
dR1_dtheta_numerical = numerical_derivative(lambda theta: compute_R1(F_test, theta), theta_test)
print(f"dR1/dtheta numerical: {dR1_dtheta_numerical:.6e}")

# Test dR2/dF
dR2_dF_numerical = numerical_derivative(lambda F: compute_R2(F, theta_test), F_test)
print(f"dR2/dF numerical: {dR2_dF_numerical:.6e}")

# Test dR2/dtheta
dR2_dtheta_numerical = numerical_derivative(lambda theta: compute_R2(F_test, theta), theta_test)
print(f"dR2/dtheta numerical: {dR2_dtheta_numerical:.6e}")

print()

# Test second derivatives
d2R1_dF2_numerical = numerical_second_derivative(lambda F: compute_R1(F, theta_test), F_test)
print(f"d2R1/dF2 numerical: {d2R1_dF2_numerical:.6e}")

d2R1_dtheta2_numerical = numerical_second_derivative(lambda theta: compute_R1(F_test, theta), theta_test)
print(f"d2R1/dtheta2 numerical: {d2R1_dtheta2_numerical:.6e}")

d2R2_dF2_numerical = numerical_second_derivative(lambda F: compute_R2(F, theta_test), F_test)
print(f"d2R2/dF2 numerical: {d2R2_dF2_numerical:.6e}")

d2R2_dtheta2_numerical = numerical_second_derivative(lambda theta: compute_R2(F_test, theta), theta_test)
print(f"d2R2/dtheta2 numerical: {d2R2_dtheta2_numerical:.6e}")

# Mixed derivative
def compute_R2_mixed(F, theta):
    return compute_R2(F, theta)

h = 1e-6
d2R2_dFdtheta_numerical = (compute_R2_mixed(F_test + h, theta_test + h) - compute_R2_mixed(F_test + h, theta_test - h) - compute_R2_mixed(F_test - h, theta_test + h) + compute_R2_mixed(F_test - h, theta_test - h)) / (4 * h**2)
print(f"d2R2/dFdtheta numerical: {d2R2_dFdtheta_numerical:.6e}")