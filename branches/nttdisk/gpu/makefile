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
	$(CC) $(CFLAGS) -DMAIN modular_arithmetic.c $(LDFLAGS) -o $@
stage1-c.o : stage1-c.c
	$(CC) $(CFLAGS) -c stage1-c.c -o stage1-c.o
mod.o: modular_arithmetic.c modular_arithmetic.h
	$(CC) $(CFLAGS) -c modular_arithmetic.c -o mod.o
stage1-c : stage1-c.o getprime.o mod.o
	$(CC) $(CFLAGS) stage1-c.o mod.o getprime.o -o stage1-c $(LDFLAGS)
clean:
	\rm -f stage1 modular_arithmetic stage1-c *.o
