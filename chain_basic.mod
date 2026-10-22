# constraints 
param m integer > 1;
param L > 0;

# init final cords
param a;
param b;

# short
param c {k in 1..m} := m - k + 0.5;

# init range
var x {1..m};
var y {1..m};

# min of sum
minimize PotentialEnergy:
    sum {k in 1..m} c[k] * y[k];

# length of each section
subject to LinkLength {k in 1..m}:
    x[k]^2 + y[k]^2 = (L/m)^2;

# midpoint pos
subject to EndPointX:
    sum {k in 1..m} x[k] = a;

subject to EndPointY:
    sum {k in 1..m} y[k] = b;

