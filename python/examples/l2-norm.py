import diskload
from math import sqrt
import random

W = 1.0

trials = 100

love = diskload.love_numbers.read()
love = diskload.love_numbers.extrapolate( love, 4000000 )
  
csv = open('l2-norm.csv', 'w')
csv.write("cutoff U V G\n")

for cutoff in range(40000,4000000,250000):
    uTotal = 0.0;
    vTotal = 0.0;
    gTotal = 0.0;    

    for _ in range(trials):
        alpha = random.uniform(0,1)
        theta = random.uniform(0,1)
        uC, vC, gC = diskload.truncated( alpha, diskload.Compensation.UNCOMPENSATED, theta, W, cutoff, love )
        u, v, g = diskload.elliptic( alpha, diskload.Compensation.UNCOMPENSATED, theta, W, 40000, love )
        uTotal += ( uC - u )*( uC - u );
        vTotal += ( vC - v )*( vC - v );
        gTotal += ( gC - g )*( gC - g );      

    line = "%d %e %e %e\n" % (cutoff, sqrt(uTotal/trials),sqrt(vTotal/trials),sqrt(gTotal/trials))
    csv.write(line)
    print(line)
    
