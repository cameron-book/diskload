set terminal tikz
set output 'l2-norm.tex'
set logscale y
set logscale x
set style line 1 dt 1 lc 'black'
set style line 2 dt 2 lc 'black'
set style line 3 dt 3 lc 'black'
set key right top
set format y "$10^{%T}$"
set title '$L^2$ norm of difference between core and truncation methods'
set xlabel 'terms in truncated series'
set ylabel '$L^2$ norm of difference'
plot 'l2-norm.csv' using (column("cutoff")):(column("U")) with lines ls 1 title "$U$", 'l2-norm.csv' using (column("cutoff")):(column("V")) with lines ls 2 title "$V$", 'l2-norm.csv' using (column("cutoff")):(column("G")) with lines ls 3 title "$G$"
