objects := $(wildcard *.o)
libs = -lm
FLAGS = -g

heat : heat.o grid.o field.o utils.o
	cc -O3 -g -o $@ heat.o grid.o field.o utils.o $(libs)


all : heat

clean : 
	rm -f $(objects) heat

%.o : %.c
	cc -O3 -g -c $(CFLAGS) $<
