CC=gcc
CFLAGS=-O2 -g -W -Wall -pedantic -ansi -I$(GMP)/include
LDFLAGS=$(GMP)/lib/libgmp.a -lm

stage1 : stage1.o prototype.o getprime.o
	$(CC) $(CFLAGS) stage1.o prototype.o getprime.o -o stage1 $(LDFLAGS)
stage1.o : stage1.c prototype.h getprime.h
	$(CC) $(CFLAGS) -c stage1.c -o stage1.o
prototype.o : prototype.c prototype.h
	$(CC) $(CFLAGS) -c prototype.c -o prototype.o
getprime.o : getprime.c getprime.h
	$(CC) $(CFLAGS) -c getprime.c -o getprime.o

modular_arithmetic: modular_arithmetic.c
	$(CC) $(CFLAGS) modular_arithmetic.c $(LDFLAGS) -o $@
clean:
	\rm -f stage1 modular_arithmetic *.o
