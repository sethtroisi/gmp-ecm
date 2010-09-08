CC=gcc
CFLAGS=-O2 -g -W -Wall -pedantic -ansi

stage1 : stage1.o prototype.o getprime.o
	$(CC) $(CFLAGS) stage1.o prototype.o getprime.o -o stage1 -lgmp -lm
stage1.o : stage1.c prototype.h getprime.h
	$(CC) $(CFLAGS) -c stage1.c -o stage1.o -lgmp -lm
prototype.o : prototype.c prototype.h
	$(CC) $(CFLAGS) -c prototype.c -o prototype.o -lgmp -lm
getprime.o : getprime.c getprime.h
	$(CC) $(CFLAGS) -c getprime.c -o getprime.o -lgmp -lm

a.out: modular_arithmetic.c
	$(CC) $(CFLAGS) modular_arithmetic.c -lgmp -lm
clean:
	\rm stage1 a.out *.o
