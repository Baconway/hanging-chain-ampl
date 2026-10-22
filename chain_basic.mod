# ==============================================================================
# Model: chain_basic.mod
# Exercise 4: The Hanging Chain Problem (Question 2 - Basic 3-Link Formulation)
# ==============================================================================

# Number of rigid links
param m integer > 1;

# Total length of the chain
param L > 0;

# Endpoints coordinates: chain goes from (0,0) to (a,b)
param a;
param b;

# Potential energy weight coefficients: c[k] = m - k + 0.5
param c {k in 1..m} := m - k + 0.5;

# Decision variables: horizontal and vertical increments along each link
var x {1..m};
var y {1..m};

# Objective function: Total gravitational potential energy
minimize PotentialEnergy:
    sum {k in 1..m} c[k] * y[k];

# Constraints:
# 1. Equal rigid link length constraint: x_k^2 + y_k^2 = (L/m)^2
subject to LinkLength {k in 1..m}:
    x[k]^2 + y[k]^2 = (L/m)^2;

# 2. Horizontal endpoint boundary constraint: sum x_k = a
subject to EndPointX:
    sum {k in 1..m} x[k] = a;

# 3. Vertical endpoint boundary constraint: sum y_k = b
subject to EndPointY:
    sum {k in 1..m} y[k] = b;

