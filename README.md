# fault-swath
Code for making a swath profile along a fault
Bisection method

1. Read in shapefile of the fault
2. Create a series of points along the fault at a defined spacing - this will be the final spacing of the aggregated data. Call this n - should be power2: n = 2^m
3. For each point, define a vector orthogonal to the line which the point is on (defined by the neighbouring points).
       --Calculate equation for line and save parameters in array (a, b, c) ax + by - c = 0 where a > 0.
4. For each point in the river dataset
       -- Start with the middle baseline point (q = 2^(m-1)).  Evaluate alpha: ax + by - c = alpha for this river point (x and y from river point, a, b, c from line)
       -- If alpha > 0, river point is to the right of the line. If alpha < 0, river point is to the left. If alpha = 0, point is on the line.
       -- Take further baseline point depending on whether alpha is positive or negative. If positive: q_new = q + 2^(m-1-i), if negative: q_new = q - 2^(m-1-i)
       -- CONTINUE until m-1-i = 0. 
       
       
