objects := $(wildcard *.o)
libs = -lm
FLAGS = -g
op = -O3

CLANG_OMP = ~/openmp/build/bin/clang

heat : heat.o grid.o field.o utils.o
	cc $(op) -g -o $@ heat.o grid.o field.o utils.o $(libs)

heat_omp : heat_omp.o grid.o field.o utils.o
	$(CLANG_OMP) -D OPENMP -fopenmp  $(op) -o $@ heat_omp.o grid.o field.o utils.o $(libs)

all : heat heat_omp

clean : 
	rm -f $(objects) heat heat_omp

heat_omp.o :
	cc $(op) -D OPENMP -o heat_omp.o -c $(CFLAGS) heat.c
%.o : %.c
	cc $(op)  -c $(CFLAGS) $<
