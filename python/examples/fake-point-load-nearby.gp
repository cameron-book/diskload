set terminal tikz
set output 'fake-point-load-nearby.tex'
set logscale x
set key right bottom
set title '$U(0.001\degree,\theta)$ versus a point load'
set style line 1 dt 1 lc 'black'
set style line 2 dt 3 lc 'black' lw 3
set style line 3 dt 3 lc 'black'
set xlabel '$\theta$ (degrees)'
set format x "$%g\\degree$"
set ylabel 'displacement (mm)'
plot 'fake-point-load-nearby.csv' using (column("theta")):(column("Upoint")) with lines title "Point load" ls 1, 'fake-point-load-nearby.csv' using (column("theta")):(column("U")) with lines title "Disk load" ls 2, 'fake-point-load-nearby.csv' using (column("theta")):(column("U40K")) with lines title "Disk load truncated at 40k" ls 3
