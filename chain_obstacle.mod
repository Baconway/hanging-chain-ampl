# ==============================================================================
# Model: chain_obstacle.mod
# Exercise 4: The Hanging Chain Problem (Question 3 - With Obstacle)
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

# Obstacle parameters for: obst(x) = h * exp( - (x - mu)^2 / sigmasquare ) - nu
param h >= 0;
param mu;
param sigmasquare > 0;
param nu;

# Number of discretization points per link
param n integer > 0 default 20;

# Decision variables: horizontal and vertical increments along each link
var x {1..m};
var y {1..m};

# Cumulative joint/node coordinates:
# Node 0 is fixed at origin (0, 0)
# Node k is (sum_{i=1}^k x_i, sum_{i=1}^k y_i)
var lx {k in 0..m} = if k = 0 then 0.0 else sum {i in 1..k} x[i];
var ly {k in 0..m} = if k = 0 then 0.0 else sum {i in 1..k} y[i];

# Coordinate discretization along each link k in 1..m:
# Line segment between node (k-1) and node k is parameterized by i/n in [0, 1]
var cx {k in 1..m, i in 0..n} = (1 - i/n) * lx[k-1] + (i/n) * lx[k];
var cy {k in 1..m, i in 0..n} = (1 - i/n) * ly[k-1] + (i/n) * ly[k];

# Objective: Minimize total gravitational potential energy
minimize PotentialEnergy:
    sum {k in 1..m} c[k] * y[k];

# Physical constraints:
# 1. Equal rigid link lengths: x_k^2 + y_k^2 = (L/m)^2
subject to LinkLength {k in 1..m}:
    x[k]^2 + y[k]^2 = (L/m)^2;

# 2. Horizontal endpoint boundary constraint: sum x_k = a
subject to EndPointX:
    sum {k in 1..m} x[k] = a;

# 3. Vertical endpoint boundary constraint: sum y_k = b
subject to EndPointY:
    sum {k in 1..m} y[k] = b;

# 4. Obstacle avoidance: entire chain profile must lie on or above the obstacle
subject to ObstacleAvoidance {k in 1..m, i in 0..n}:
    cy[k,i] >= h * exp( - (cx[k,i] - mu)^2 / sigmasquare ) - nu;

