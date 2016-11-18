objects := $(wildcard *.o)
libs = -lm
FLAGS = -O3
#FLAGS = -g -O0 -Wall

CC_OMP = ~/openmp/build/bin/clang
CC_MPI = mpicc

heat : heat.o grid.o field.o utils.o
	cc $(FLAGS) -g -o $@ heat.o grid.o field.o utils.o $(libs)

heat_omp : heat_omp.o grid.o field.o utils.o
	$(CC_OMP) -D OPENMP -fopenmp  $(FLAGS) -o $@ heat_omp.o grid.o field.o utils.o $(libs)

heat_mpi : heat_mpi.o grid.o field.o utils.o
	$(CC_MPI) -D MPI  $(FLAGS) -o $@ heat_mpi.o grid.o field.o utils.o $(libs)

all : heat heat_omp heat_mpi

clean : 
	rm -f $(objects) heat heat_omp heat_mpi

heat_omp.o : heat.c
	$(CC_OMP) $(FLAGS) -D OPENMP -fopenmp -o heat_omp.o -c $(CFLAGS) heat.c

heat_mpi.o : heat.c
	$(CC_MPI) $(FLAGS) -D MPI -o heat_mpi.o -c $(CFLAGS) heat.c
%.o : %.c
	cc $(FLAGS) -c $(CFLAGS) $<
