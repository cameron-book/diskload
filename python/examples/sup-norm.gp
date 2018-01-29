set terminal tikz
set output 'sup-norm.tex'
set logscale y
set style line 1 dt 1 lc 'black'
set style line 2 dt 2 lc 'black'
set style line 3 dt 3 lc 'black'
set key right top
set format y "$10^{%T}$"
set title 'Difference between core and truncation away from disk edge'
set xlabel 'terms in truncated series'
set ylabel 'sup norm of difference'
plot 'sup-norm.csv' using (column("cutoff")):(column("U")) with lines ls 1 title "$U$", 'sup-norm.csv' using (column("cutoff")):(column("V")) with lines ls 2 title "$V$", 'sup-norm.csv' using (column("cutoff")):(column("G")) with lines ls 3 title "$G$"
