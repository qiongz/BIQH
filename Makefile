CC = icpc
CFLAGS =-std=c++11 -DMKL_ILP64 -mkl=parallel -qopenmp
LIBS=   -liomp5 -lpthread -lm -ldl

blqh:main.cpp basis.o matrix.o init.o hamiltonian.o mt19937-64.o lanczos_hamiltonian.o
	$(CC) $(CFLAGS) $^ -O3 -o $@ ${LIBS} $(CFLAGS) -lgsl

basis.o:basis.cpp basis.h
	$(CC) $(CFLAGS) -c basis.cpp

matrix.o:matrix.cpp matrix.h mt19937-64.h
	$(CC) $(CFLAGS) -c matrix.cpp -o $@

init.o:init.cpp init.h
	$(CC) $(CFLAGS) -c init.cpp

hamiltonian.o:hamiltonian.cpp hamiltonian.h matrix.h
	$(CC) $(CFLAGS) -c hamiltonian.cpp

lanczos_hamiltonian.o:lanczos_hamiltonian.cpp lanczos_hamiltonian.h matrix.h
	$(CC) $(CFLAGS) -c lanczos_hamiltonian.cpp

mt19937-64.o:mt19937-64.c mt19937-64.h
	$(CC) -c mt19937-64.c

.PHONY: all blqh remove
all: clean blqh 

clean:
	rm -f  *.o
remove:
	rm  blqh
