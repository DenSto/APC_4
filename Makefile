objects := $(wildcard *.o)
libs = -lm
FLAGS = -g

heat : heat.o grid.o field.o
	cc -O3 -o $@ heat.o grid.o field.o $(libs)


all : heat

clean : 
	rm -f $(objects) heat

%.o : %.c
	cc -O3 -c $(CFLAGS) $<
