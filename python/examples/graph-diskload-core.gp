set terminal tikz
set output 'graph-diskload-core.tex'
set key right bottom
set title 'Disk load computed via core'
set style line 1 dt 1 lc 'black'
set style line 2 dt 2 lc 'black'
set style line 3 dt 3 lc 'black'
set xlabel '$\theta/\alpha$'
set ylabel 'displacement (mm)'
plot 'wide-data.csv' using (column("alpha/theta")):(column("U")) with lines ls 1 title "$U$",  'wide-data.csv' using (column("alpha/theta")):(column("V")) with lines ls 2 title "$V$",  'wide-data.csv' using (column("alpha/theta")):(column("G")) with lines ls 3 title "$G$"
