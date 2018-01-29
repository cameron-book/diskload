import diskload
import math
from math import pi, sin
import random

uncompensated = diskload.Compensation.UNCOMPENSATED
love = diskload.love_numbers.read()
love = diskload.love_numbers.extrapolate( love, 4000000 )

N = 1000

alpha = 0.001
min_theta = alpha * 10.0
max_theta = alpha * 1000.0

base = math.exp( math.log( max_theta / min_theta ) / (N-1) )

f = open("fake-point-load.csv", "w")
f.write( "alpha theta alpha/theta U V G U40K V40K G40K Upoint Vpoint Gpoint\n" )

for i in range(N):
    theta = (base ** i) * min_theta

    w = 1.0
    
    f.write("%e %e %e " % (alpha, theta, theta/alpha) )
    f.write("%e %e %e " % diskload.elliptic( alpha, uncompensated, theta, w, 40000, love ) )
    f.write("%e %e %e " % diskload.truncated( alpha, uncompensated, theta, w, 40000, love ) )

    radius = 2.0 * 1000.0 * diskload.earth.default.radius * sin(alpha*diskload.constants.DEGREES / 2.0)
    area = (pi * radius * radius)
    w = area * 4.0 * diskload.earth.default.gravity / 1e16
    
    f.write("%e %e %e " % diskload.point( theta, w, 40000, love ) )
    f.write("\n")
