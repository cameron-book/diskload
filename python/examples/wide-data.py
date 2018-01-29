import diskload
import math
import random

uncompensated = diskload.Compensation.UNCOMPENSATED
love = diskload.love_numbers.read()
love = diskload.love_numbers.extrapolate( love, 4000000 )

N = 100
w = 1.0

alpha = 0.1
min_theta = alpha * 0.5
max_theta = alpha * 2.0

base = math.exp( math.log( max_theta / min_theta ) / (N-1) )

f = open("wide-data.csv", "w")

f.write("alpha theta alpha/theta U4K V4K G4K U40K V40K G40K U V G U400K V400K G400K U4M V4M G4M\n" )
  
for i in range(N):
    theta = (base ** i) * min_theta

    f.write("%e %e %e " % (alpha, theta, theta/alpha) )
    
    f.write("%e %e %e " % diskload.truncated( alpha, uncompensated, theta, w, 4000, love ) )
    f.write("%e %e %e " % diskload.truncated( alpha, uncompensated, theta, w, 40000, love ) )
    f.write("%e %e %e " % diskload.elliptic( alpha, uncompensated, theta, w, 40000, love ) )
    f.write("%e %e %e " % diskload.truncated( alpha, uncompensated, theta, w, 400000, love ) )
    f.write("%e %e %e " % diskload.truncated( alpha, uncompensated, theta, w, 4000000, love ) )    
    f.write("\n")
