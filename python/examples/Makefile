l2-norm.csv: l2-norm.py
	PYTHONPATH=.. python l2-norm.py

data.csv: data.py
	PYTHONPATH=.. python data.py

wide-data.csv: wide-data.py
	PYTHONPATH=.. python wide-data.py

l2-norm.tex: l2-norm.gp l2-norm.csv
	gnuplot l2-norm.gp

timings.tex table-timings.tex timings.gp timings.csv: timings.py
	PYTHONPATH=.. python timings.py
	gnuplot timings.gp

timings-u.tex table-timings-u.tex timings-u.gp timings-u.csv: timings.py
	PYTHONPATH=.. python timings.py u
	gnuplot timings-u.gp

timings-v.tex table-timings-v.tex timings-v.gp timings-v.csv: timings.py
	PYTHONPATH=.. python timings.py v
	gnuplot timings-v.gp

timings-g.tex table-timings-g.tex timings-g.gp timings-g.csv: timings.py
	PYTHONPATH=.. python timings.py g
	gnuplot timings-g.gp

fake-point-load.csv: fake-point-load.py
	PYTHONPATH=.. python fake-point-load.py 

fake-point-load-difference.tex: fake-point-load-difference.gp fake-point-load.csv
	gnuplot fake-point-load-difference.gp

fake-point-load.tex: fake-point-load.gp fake-point-load.csv
	gnuplot fake-point-load.gp 

%.tex: %.gp data.csv wide-data.csv
	gnuplot $<

assets.pdf: assets.tex graph-core-versus-400k.tex graph-core-versus-40k.tex graph-core-versus-4k.tex graph-core-versus-4m.tex graph-diskload-core.tex l2-norm.tex table-timings-u.tex table-timings-v.tex table-timings-g.tex timings-u.tex timings-v.tex timings-g.tex table-timings.tex timings.tex fake-point-load-difference.tex fake-point-load.tex
	pdflatex assets.tex

clean:
	rm -f *.aux *.pdf *.csv *.log timings-g.tex table-timings-g.tex timings-g.gp timings-g.csv timings-v.tex table-timings-v.tex timings-v.gp timings-v.csv timings-u.tex table-timings-u.tex timings-u.gp timings-u.csv timings.tex table-timings.tex timings.gp timings.csv 
	rm -f graph-core-versus-400k.tex graph-core-versus-40k.tex graph-core-versus-4k.tex graph-core-versus-4m.tex graph-diskload-core.tex l2-norm.tex 


