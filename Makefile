objects := $(wildcard *.o)
libs = -lm
FLAGS = -g
op = -O3
CLANG_OMP = ~/openmp/build/bin/clang

heat : heat.o grid.o field.o utils.o
	cc $(op) -g -o $@ heat.o grid.o field.o utils.o $(libs)
heat_omp : heat_omp.o grid.o field.o utils.o
	$(CLANG_OMP) -fopenmp  $(op) -o $@ heat_omp.o grid.o field.o utils.o $(libs)


all : heat heat_omp

clean : 
	rm -f $(objects) heat heat_omp

%.o : %.c
	$(CLANG_OMP) -fopenmp $(op) -c $(CFLAGS) $<
#%.o : %.c
#	cc -O3 -g -c $(CFLAGS) $<
