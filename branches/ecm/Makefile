GMPDIR=/global/ecrouves/loria/linux/gmp-4.1
I=$(GMPDIR)/include
L=$(GMPDIR)/lib
CFLAGS=-g -Wall -Wmissing-prototypes -ansi -pedantic
CC=gcc

FILES= aux.o bestd.o ecm.o getprime.o main.o pm1.o poly.o stage2.o
DIST= COPYING Makefile README aux.c bestd.c check.mpl cputime.h ecm.c ecm.h getprime.c main.c pm1.c poly.c stage2.c

.SUFFIXES: .c .o

ecm: $(FILES)
	$(CC) $(CFLAGS) -L$(L) $(FILES) -o $@ -lgmp -lm

.c.o:
	$(CC) $(CFLAGS) -I$(I) -c $<

clean:
	rm ecm *.o *~

dist: $(DIST)
	mkdir ecm-5.0
	cp $(DIST) ecm-5.0
	tar cf ecm-5.0.tar ecm-5.0
	/bin/rm -fr ecm-5.0
	gzip --best ecm-5.0.tar
