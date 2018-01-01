elliptic-integrals-test: elliptic-integrals-test.c elliptic-integrals.c tanhsinh.c
	gcc -O2 -o elliptic-integrals-test elliptic-integrals-test.c elliptic-integrals.c tanhsinh.c `pkg-config --cflags --libs gsl`
	./elliptic-integrals-test

test: test.c diskload.c elliptic-integrals.c tanhsinh.c
	gcc -O2 -o test diskload.c test.c elliptic-integrals.c tanhsinh.c `pkg-config --cflags --libs gsl`
	./test

benchmark: benchmark.c diskload.c elliptic-integrals.c tanhsinh.c
	gcc -no-pie -pg -o benchmark diskload.c elliptic-integrals.c tanhsinh.c benchmark.c `pkg-config --cflags --libs gsl`
	#gcc -O2 -o benchmark diskload.c elliptic-integrals.c tanhsinh.c benchmark.c `pkg-config --cflags --libs gsl`
	./benchmark

integrands: integrands.c diskload.c elliptic-integrals.c tanhsinh.c
	gcc -O2 -o integrands diskload.c elliptic-integrals.c tanhsinh.c integrands.c `pkg-config --cflags --libs gsl`
	./integrands

fake-point-load-difference.tex: fake-point-load-difference.gp fake-point-load.csv
	gnuplot fake-point-load-difference.gp 

fake-point-load.csv: fake-point-load.c diskload.c elliptic-integrals.c tanhsinh.c
	gcc -O2 -o fake-point-load diskload.c elliptic-integrals.c tanhsinh.c fake-point-load.c `pkg-config --cflags --libs gsl`
	./fake-point-load

fake-point-load-nearby.csv: fake-point-load.c diskload.c elliptic-integrals.c tanhsinh.c
	gcc -DPOSTFIX=\"-nearby\" -DMIN_THETA=1.0 -DMAX_THETA=10.0 -O2 -o fake-point-load diskload.c elliptic-integrals.c tanhsinh.c fake-point-load.c `pkg-config --cflags --libs gsl`
	./fake-point-load

timings.csv timings.gp table-timings.tex: timings.c diskload.c elliptic-integrals.c tanhsinh.c
	gcc -O2 -o timings diskload.c elliptic-integrals.c tanhsinh.c timings.c `pkg-config --cflags --libs gsl`
	./timings

timings.tex: timings.gp timings.csv
	gnuplot timings.gp

timings-u.csv timings-u.gp table-timings-u.tex: timings.c diskload.c elliptic-integrals.c tanhsinh.c
	gcc -DPOSTFIX=\"-u\" -DTITLE='" for $$U$$"' -DDISKLOAD_NO_V -DDISKLOAD_NO_G -O2 -o timings diskload.c elliptic-integrals.c tanhsinh.c timings.c `pkg-config --cflags --libs gsl`
	./timings

timings-u.tex: timings-u.gp timings-u.csv
	gnuplot timings-u.gp

timings-v.csv timings-v.gp table-timings-v.tex: timings.c diskload.c elliptic-integrals.c tanhsinh.c
	gcc -DPOSTFIX=\"-v\" -DTITLE='" for $$V$$"' -DDISKLOAD_NO_U -DDISKLOAD_NO_G -O2 -o timings diskload.c elliptic-integrals.c tanhsinh.c timings.c `pkg-config --cflags --libs gsl`
	./timings

timings-v.tex: timings-v.gp timings-v.csv
	gnuplot timings-v.gp

timings-g.csv timings-g.gp table-timings-g.tex: timings.c diskload.c elliptic-integrals.c tanhsinh.c
	gcc -DPOSTFIX=\"-g\" -DTITLE='" for $$G$$"' -DDISKLOAD_NO_U -DDISKLOAD_NO_V -O2 -o timings diskload.c elliptic-integrals.c tanhsinh.c timings.c `pkg-config --cflags --libs gsl`
	./timings

timings-g.tex: timings-g.gp timings-g.csv
	gnuplot timings-g.gp


sup-norm.tex: sup-norm.gp sup-norm.csv
	gnuplot sup-norm.gp

l2-norm.tex: l2-norm.gp l2-norm.csv
	gnuplot l2-norm.gp

%.tex: %.gp data.csv
	gnuplot $<

data.csv: data.c diskload.c elliptic-integrals.c tanhsinh.c
	gcc -o data diskload.c elliptic-integrals.c tanhsinh.c data.c `pkg-config --cflags --libs gsl`
	./data

sup-norm.csv: sup-norm.c diskload.c elliptic-integrals.c tanhsinh.c
	gcc -o sup-norm diskload.c elliptic-integrals.c tanhsinh.c sup-norm.c `pkg-config --cflags --libs gsl`
	./sup-norm

l2-norm.csv: l2-norm.c diskload.c elliptic-integrals.c tanhsinh.c
	gcc -o l2-norm diskload.c elliptic-integrals.c tanhsinh.c l2-norm.c `pkg-config --cflags --libs gsl`
	./l2-norm

assets.pdf: assets.tex sup-norm.tex timings.tex table-timings.tex graph-core-versus-400k.tex graph-core-versus-40k.tex graph-core-versus-4k.tex graph-core-versus-4m.tex  graph-diskload-core.tex
	pdflatex assets.tex
