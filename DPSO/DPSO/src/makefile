# Makefile for GNU make

LDLIBS = -lm 
CC = gcc
CFLAGS =  -O2  -Wall

.SUFFIXES: .o .cc .c

OBJS_c = pso.c myfun.c methods.c mem_loc.c
OBJS_o = pso.o myfun.o methods.o mem_loc.o

dpso:   $(OBJS_c) $(OBJS_o)
	$(CC) $(LDLIBS) $(CFLAGS) -o $@ $(OBJS_o) $(LDLIBS)

