objects := $(wildcard *.o)
libs = -lm
FLAGS = -O3
#FLAGS = -g -O0 -Wall

CLANG_OMP = ~/openmp/build/bin/clang
CLANG_MPI = mpicc

heat : heat.o grid.o field.o utils.o
	cc $(FLAGS) -g -o $@ heat.o grid.o field.o utils.o $(libs)

heat_omp : heat_omp.o grid.o field.o utils.o
	$(CLANG_OMP) -D OPENMP -fopenmp  $(FLAGS) -o $@ heat_omp.o grid.o field.o utils.o $(libs)

heat_mpi : heat.c heat_mpi.o grid.o field.o utils.o
	$(CLANG_MPI) -D MPI  $(FLAGS) -o $@ heat_mpi.o grid.o field.o utils.o $(libs)

all : heat heat_omp heat_mpi

clean : 
	rm -f $(objects) heat heat_omp heat_mpi

heat_omp.o :
	cc $(FLAGS) -D OPENMP -o heat_omp.o -c $(CFLAGS) heat.c
heat_mpi.o :
	$(CLANG_MPI) $(FLAGS) -D MPI -o heat_mpi.o -c $(CFLAGS) heat.c
%.o : %.c
	cc $(FLAGS) -c $(CFLAGS) $<
