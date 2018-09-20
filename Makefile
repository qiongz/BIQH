CC               := g++
ICC              := icpc
CFLAGS           := -std=c++11
ICC_FLAGS        := -std=c++11 -qopenmp -Dmkl
LDFLAGS          := -O2 -fopenmp -lpthread
ICC_LDFLAGS      := -O2 -qopenmp
LAPACK_INC       :=
LAPACK_LIBS      := -llapack -lblas
MKL_INC          := -m64 -I$(MKLROOT)/include 
SEQ_MKL_LIBS     := -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 \
       	  	     -lmkl_sequential -lmkl_core -lpthread -lm -ldl 

PARAL_MKL_LIBS   := -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 \
	             -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl


all:blqh blqh_mkl clean

.PHONY:clean remove
clean:
	rm -f *.o 

remove:
	rm -f blqh blqh_mkl

#------------- lapack version ---------------#

blqh:main.o  basis.o matrix.o init.o  mt19937-64.o lanczos_hamiltonian.o hamiltonian.o
	$(CC) $(LDFLAGS) $^ -o $@  $(LAPACK_LIBS)  -lgsl

main.o:main.cpp
	$(CC) $(CFLAGS) $(LAPACK_INC) -c -o $@ $<

basis.o:basis.cpp  basis.h
	$(CC) $(CFLAGS) -c -o $@ $<

matrix.o:matrix.cpp matrix.h mt19937-64.h
	$(CC) $(CFLAGS) $(LAPACK_INC) -c -o $@ $<

init.o:init.cpp init.h
	$(CC) $(CFLAGS) -c -o $@ $<

hamiltonian.o:hamiltonian.cpp hamiltonian.h matrix.h
	$(CC) $(CFLAGS) $(LAPACK_INC) -c -o $@ $<

lanczos_hamiltonian.o:lanczos_hamiltonian.cpp lanczos_hamiltonian.h matrix.h
	$(CC) $(CFLAGS) $(LAPACK_INC) -c -o $@ $<

mt19937-64.o:mt19937-64.c mt19937-64.h
	$(CC) -c -o $@ $<

#----------------- ICC MKL version ----------#

blqh_mkl:main_icc.o  basis_icc.o matrix_icc.o init_icc.o lanczos_hamiltonian_icc.o hamiltonian_icc.o
	$(ICC) $(ICC_LDFLAGS) $^ -o $@  $(PARAL_MKL_LIBS)  -lgsl

main_icc.o:main.cpp
	$(ICC) $(ICC_FLAGS) $(MKL_INC) -c -o $@ $<

basis_icc.o:basis.cpp  basis.h
	$(ICC) $(ICC_FLAGS) -c -o $@ $<

matrix_icc.o:matrix.cpp matrix.h mt19937-64.h
	$(ICC) $(ICC_FLAGS) $(MKL_INC) -c -o $@ $<

init_icc.o:init.cpp init.h
	$(ICC) $(ICC_FLAGS) -c -o $@ $<

hamiltonian_icc.o:hamiltonian.cpp hamiltonian.h matrix.h
	$(ICC) $(ICC_FLAGS) $(MKL_INC) -c -o $@ $<

lanczos_hamiltonian_icc.o:lanczos_hamiltonian.cpp lanczos_hamiltonian.h matrix.h
	$(ICC) $(ICC_FLAGS) $(MKL_INC) -c -o $@ $<
