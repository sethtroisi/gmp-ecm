# Makefile for gmp-ecm
#
# Copyright 2001, 2002, 2003 Alexander Kruppa and Paul Zimmermann.
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# directory where GMP is installed
# gmp.h should be in $(GMP)/include
# libgmp.a/libgmp.so should be in $(GMP)/lib
GMP=/usr/local/gmp

# directory where NTL is installed
# ZZ_pX.h and version.h should be in $(NTL)/include/NTL
# libntl.a should be in $(NTL)/lib
NTL=/usr/local/ntl

VERSION=5.0-beta2

###################### do not edit below this line ############################

FILES= auxi.o bestd.o ecm.o ecm2.o getprime.o listz.o lucas.o main.o pm1.o pp1.o stage2.o toomcook.o memory.o mpmod.o mul_lo.o polyeval.o resume.o

CFLAG= -O2 -g -Wall -Wmissing-prototypes
LDFLAG= -lgmp -lm
CXX=gcc
CC=gcc
LD=gcc
EXTRAFILES=
ALLFILES= $(FILES) $(EXTRAFILES)
CFLAGS= $(CFLAG)
LDFLAGS= $(LDFLAG)
POLYGCD=0

DIST=  auxi.c bestd.c ecm.c ecm2.c getprime.c listz.c lucas.c main.c ntl.c pm1.c polyz.c pp1.c stage2.c toomcook.c memory.c mpmod.c mul_lo.c polyeval.c resume.c
EXTRADIST= COPYING INSTALL Makefile README ecm.h test.pm1 test.pp1 test.ecm tune.c

.SUFFIXES: .c .o

all:
	@if test $(POLYGCD) -ne 0; then \
           CXX=g++;                    \
           EXTRAFILES=ntl.o;           \
           LDFLAGS="-lntl $(LDFLAG)";  \
           make ecm GMP=$(GMP) NTL=$(NTL) CXX=g++ EXTRAFILES="ntl.o polyz.o" LDFLAGS="-lntl $(LDFLAG) CFLAGS="$(CFLAG) -DPOLYGCD""; \
        else \
           make ecm GMP=$(GMP); \
        fi

ecm: $(ALLFILES) ecm.h
	$(LD) $(CFLAGS) -L$(GMP)/lib -L$(NTL)/lib $(ALLFILES) -o $@ $(LDFLAGS)

tune: mpmod.o ecm.h tune.o auxi.o mul_lo.o
	$(CC) $(CFLAGS) -L$(GMP)/lib tune.o mpmod.o auxi.o mul_lo.o -o $@ -lgmp

ntl.o: ntl.c
	$(CXX) $(CFLAGS) -c -I$(GMP)/include -I$(NTL)/include $<

.c.o: ecm.h
	$(CC) $(CFLAGS) -I$(GMP)/include -c $<

clean:
	rm ecm *.o *~

dist: $(DIST)
	mkdir ecm-$(VERSION)
	cp $(DIST) $(EXTRADIST) ecm-$(VERSION)
	tar cf ecm-$(VERSION).tar ecm-$(VERSION)
	/bin/rm -fr ecm-$(VERSION)
	gzip --best ecm-$(VERSION).tar
