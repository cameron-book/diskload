set terminal tikz
set output 'fake-point-load-difference.tex'
set logscale x
set logscale y
set format y "$10^{%T}$"
unset key
set title '$U(0.001\degree,\theta)$ versus a point load'
set style line 1 dt 1 lc 'black'
set style line 2 dt 3 lc 'black' lw 3
set style line 3 dt 3 lc 'black'
set xlabel '$\theta$ (degrees)'
set format x "$%g\\degree$"
set ylabel 'difference in displacement (mm)'
plot 'fake-point-load.csv' using (column("theta")):(abs(column("U") - column("Upoint"))) with lines ls 1
