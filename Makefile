# directory where GMP is installed (include, lib)
GMPDIR=/global/ecrouves/loria/linux/gmp-4.1

# directory where NTL is installed (include, lib)
NTL=/global/ecrouves/loria/linux/ntl-5.3

###################### do not edit below this line ############################

CFLAGS=-O2 -g -Wall -Wmissing-prototypes -ansi -pedantic
CC=g++

FILES= aux.o bestd.o ecm.o getprime.o main.o pm1.o listz.o stage2.o polyz.o ntl.o toomcook.o
DIST= COPYING Makefile README aux.c bestd.c check.mpl cputime.h ecm.c ecm.h getprime.c main.c pm1.c listz.c stage2.c polyz.c ntl.c toomcook.c

.SUFFIXES: .c .o

ecm: $(FILES) ecm.h
	$(CC) $(CFLAGS) -L$(GMPDIR)/lib -L$(NTL)/lib $(FILES) -o $@ -lntl -lgmp -lm

ntl.o: ntl.c
	$(CC) $(CFLAGS) -c -I$(GMPDIR)/include -I$(NTL)/include ntl.c

.c.o: ecm.h
	$(CC) $(CFLAGS) -I$(GMPDIR)/include -c $<

clean:
	rm ecm *.o *~

dist: $(DIST)
	mkdir ecm-5.0
	cp $(DIST) ecm-5.0
	tar cf ecm-5.0.tar ecm-5.0
	/bin/rm -fr ecm-5.0
	gzip --best ecm-5.0.tar
