set terminal tikz
set output 'graph-core-versus-4k.tex'
set logscale x
set arrow from 1, graph 0 to 1, graph 1 nohead lc rgb 'gray'
set xtics (0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4)
set style line 1 dt 1 lc 'black'
set style line 2 dt 2 lc 'black'
set style line 3 dt 3 lc 'black'
set key right bottom
set title 'Difference between core method and truncation to 4k terms'
set xlabel '$\theta/\alpha$'
set ylabel 'difference (mm)'
plot 'data.csv' using (column("alpha/theta")):(column("U") - column("U4K")) with lines ls 1 title "$U$", 'data.csv' using (column("alpha/theta")):(column("V") - column("V4K")) with lines ls 2 title "$V$", 'data.csv' using (column("alpha/theta")):(column("G") - column("G4K")) with lines ls 3 title "$G$"
