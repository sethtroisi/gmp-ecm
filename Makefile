# directory where GMP is installed (include, lib)
GMP=/global/ecrouves/loria/linux/gmp-4.1

# directory where NTL is installed (include, lib)
NTL=/global/ecrouves/loria/linux/ntl-5.3

###################### do not edit below this line ############################

CFLAGS=-O2 -g -Wall -Wmissing-prototypes -ansi -pedantic
CC=g++

FILES= aux.o bestd.o ecm.o getprime.o listz.o main.o ntl.o pm1.o polyz.o pp1.o stage2.o toomcook.o
DIST=  aux.c bestd.c ecm.c getprime.c listz.c main.c ntl.c pm1.c polyz.c pp1.c stage2.c toomcook.c
EXTRADIST= COPYING Makefile README cputime.h ecm.h test.pm1 test.pp1

.SUFFIXES: .c .o

ecm5: $(FILES) ecm.h
	$(CC) -static $(CFLAGS) -L$(GMP)/lib -L$(NTL)/lib $(FILES) -o $@ -lntl -lgmp -lm

ntl.o: ntl.c
	$(CC) $(CFLAGS) -c -I$(GMP)/include -I$(NTL)/include ntl.c

.c.o: ecm.h
	$(CC) $(CFLAGS) -I$(GMP)/include -c $<

clean:
	rm ecm *.o *~

dist: $(DIST)
	mkdir ecm-5.0
	cp $(DIST) $(EXTRADIST) ecm-5.0
	tar cf ecm-5.0.tar ecm-5.0
	/bin/rm -fr ecm-5.0
	gzip --best ecm-5.0.tar
