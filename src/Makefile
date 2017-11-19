test: test.c diskload.c
	gcc -o test diskload.c test.c `pkg-config --cflags --libs gsl`
	./test
