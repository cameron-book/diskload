test: test.c diskload.c
	gcc -o test diskload.c test.c `pkg-config --cflags --libs gsl`
	./test

benchmark: benchmark.c diskload.c
	gcc -no-pie -pg -o benchmark diskload.c benchmark.c `pkg-config --cflags --libs gsl`
	./benchmark
