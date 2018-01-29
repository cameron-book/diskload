import diskload
from math import sqrt
import random

W = 1.0

trials = 1500

love = diskload.love_numbers.read()
love = diskload.love_numbers.extrapolate( love, 4000000 )
  
csv = open('sup-norm.csv', 'w')
csv.write("cutoff U V G\n")

for cutoff in range(20000,4000000,100000):
    uMax = 0
    vMax = 0
    gMax = 0

    for _ in range(trials):
        alpha = random.uniform(0,1)
        theta = random.uniform(0,1)
        while abs(alpha-theta) < 0.1 or alpha < 0.1 or theta < 0.1:
            alpha = random.uniform(0,1)
            theta = random.uniform(0,1)
            
        uC, vC, gC = diskload.truncated( alpha, diskload.Compensation.UNCOMPENSATED, theta, W, cutoff, love )
        u, v, g = diskload.elliptic( alpha, diskload.Compensation.UNCOMPENSATED, theta, W, 40000, love )
        if abs(uC - u) > uMax:
            uMax = abs(uC - u)
        if abs(vC - v) > vMax:
            vMax = abs(vC - v)
        if abs(gC - g) > gMax:
            gMax = abs(gC - g)            

    line = "%d %e %e %e\n" % (cutoff, uMax, vMax, gMax)
    csv.write(line)
    print(line)
    
