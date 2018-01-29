import diskload
import math
import random
import sys

love = diskload.love_numbers.read()
love = diskload.love_numbers.extrapolate( love, 4000000 )

postfix = ""
computeU = True
computeV = True
computeG = True
title = ""

if len(sys.argv) > 1:
    postfix = "-" + sys.argv[1]
    computeU = (postfix == "-u")
    computeV = (postfix == "-v")
    computeG = (postfix == "-g")    
    title = "for $" + sys.argv[1].upper() + "$"

table = open("table-timings" + postfix + ".tex","w")
interval = 0.2

table.write( "\\begin{table}[ht]\n" )
table.write( "\\centering\n" )
table.write( "\\begin{tabular}{lS[table-format=3.2]}\n" )
table.write( "method & {time (ms)} \\\\\n" )
table.write( "\\hline\n" )

timings = open("timings" + postfix + ".csv", "w")
timings.write( "cutoff milliseconds\n" )

trials = 2000
W = 1.0

def truncated(cutoff):
    alpha = random.uniform(0,interval)
    theta = random.uniform(0,interval)
    u, v, g = diskload.truncated( alpha, diskload.Compensation.UNCOMPENSATED, theta, W, cutoff, love, computeU=computeU, computeV=computeV, computeG=computeG )

def elliptic():
    alpha = random.uniform(0,interval)
    theta = random.uniform(0,interval)
    u, v, g = diskload.elliptic( alpha, diskload.Compensation.UNCOMPENSATED, theta, W, 40000, love, computeU=computeU, computeV=computeV, computeG=computeG )
    
for cutoff in range(100000, 4000000, 100000):
    import timeit
    time_spent = timeit.timeit("truncated(" + str(cutoff) + ")", setup="from __main__ import truncated, computeU, computeV, computeG",number=trials)*1000/trials
    
    timings.write( "%d %e\n" % ( cutoff, time_spent ) )

    name = str(cutoff)
    if name[-3:] == "000":
        name = name[:-3] + "\\mathrm{k}"

    print( "%.2f ms/%d-truncated diskload\n" % (time_spent,cutoff) );
    if cutoff % 500000 == 0:
        table.write( "truncated at $N = %s$ & %.2f \\\\\n" % ( name, time_spent ) )

core_time = timeit.timeit("elliptic()", setup="from __main__ import elliptic, computeU, computeV, computeG",number=trials)*1000/trials
table.write( "core method & %.2f  \\\\\n" % core_time )
print( "%.2f ms/hypergeometric diskload\n" % core_time )
    
table.write(  "\\end{tabular}\n" )
table.write(  "\\caption{Time per disk load computation " + title + " for various methods; these are computed by sampling %d random points in the range $0 \\leq \\alpha, \\theta \\leq %.2f$}\n" % (trials, interval) )
table.write(  "\\label{table:timings" + postfix + "}\n" )
table.write(  "\\end{table}\n" )

gp = open("timings" + postfix + ".gp", "w" )

gp.write( ("set terminal tikz\n"  +
            "set output 'timings" + postfix +  ".tex'\n" +
            "set xlabel 'terms in truncated series'\n" +
            "set ylabel 'time (milliseconds)'\n" +
            "set key off\n" +
            "set style line 1 dt 1 lc 'black'\n" +
            "set arrow from graph 0, first %f to graph 1, first %f nohead lc rgb 'gray'\n" +
            "set style textbox opaque noborder\n" +
            "set title 'Speed of core method and series truncation " + title + "'\n" +
            "set label 'core method' at graph 0.7, first %f boxed front\n" +
            "plot 'timings" + postfix +  ".csv' using (column('cutoff')):(column('milliseconds')) with lines ls 1 title 'series truncation'\n") % ( core_time, core_time, core_time ) )
