CFLAGS=-Wall -O3 -march=native
LDFLAGS=

iterative_solver: iterative_solver.o
	gcc-11 -fopenmp -o iterative_solver iterative_solver.o $(LDFLAGS)

iterative_solver.o: iterative_solver.c
	gcc-11 -fopenmp $(CFLAGS) -c iterative_solver.c

clean:
	rm -f ./iterative_solver *.o