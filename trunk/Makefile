# Makefile for gmp-ecm
#
# Copyright 2001, 2002, 2003 Paul Zimmermann and Alexander Kruppa.
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

# Standard installation prefix
prefix=/usr/local

# directory where GMP is installed
# gmp.h should be in $(GMP)/include
# libgmp.a/libgmp.so should be in $(GMP)/lib
GMP=$(prefix)

# directory where NTL is installed
# ZZ_pX.h and version.h should be in $(NTL)/include/NTL
# libntl.a should be in $(NTL)/lib
NTL=$(prefix)
# example of TUNEFLAGS for x86
# TUNEFLAGS= -O3 -mcpu=i586 -fstrict-aliasing -I/usr/include/gmp -DWANT_GMP_IMPL -DMPZMOD_THRESHOLD=69 -DREDC_THRESHOLD=92
TUNEFLAGS= -O2

VERSION=5.1.2-beta

###################### do not edit below this line ############################

OBJS= auxi.o b1_ainc.o bestd.o candi.o ecm.o ecm2.o eval.o getprime.o listz.o lucas.o main.o pm1.o pp1.o stage2.o toomcook.o trial.o memory.o mpmod.o mul_lo.o polyeval.o resume.o median.o smartprp.o
CFLAGS= -g -W -Wall -Wmissing-prototypes -pedantic $(TUNEFLAGS)
LDFLAGS= -lm
CXX=g++
CC=gcc
LD=$(CC)
EXTRAOBJS=
ALLOBJS= $(OBJS) $(EXTRAOBJS)
POLYGCD=0

DIST=  auxi.c b1_ainc.c bestd.c bestdaux.c candi.c countsmooth.c ecm.c ecm2.c eval.c getprime.c listz.c lucas.c main.c median.c median-aux.c memory.c mpmod.c mul_lo.c ntl.c pm1.c polyeval.c polyz.c pp1.c resume.c smartprp.c stage2.c trial.c toomcook.c tune.c
EXTRADIST= COPYING COPYING.LIB INSTALL Makefile README ecm.h test.pm1 test.pp1 test.ecm tune.c c155 ecm-gmp.h ChangeLog

.SUFFIXES: .c .o

ecm: $(ALLOBJS) ecm.h ecm-gmp.h
	$(LD) $(ALLOBJS) -lgmp -static -o $@ $(LDFLAGS)

ecm_with_ntl:
	make ecm GMP=$(GMP) NTL=$(NTL) LD='$(CXX)' EXTRAOBJS="ntl.o polyz.o" LDFLAGS='-L$(NTL)/lib -lntl $(LDFLAGS)' CFLAGS='$(CFLAGS) -DPOLYGCD'

tune: mpmod.o ecm.h tune.o auxi.o mul_lo.o ecm-gmp.h
	$(CC) $(CFLAGS) -L$(GMP)/lib tune.o mpmod.o auxi.o mul_lo.o -o $@ $(LDFLAGS)

countsmooth: getprime.o countsmooth.o
	$(CC) $(CFLAGS) $^ -o $@ -lm -lgmp -static

bestdaux: bestdaux.c
	$(CC) $(CFLAGS) $< -o $@

ntl.o: ntl.c
	$(CXX) $(CFLAGS) -c -I$(GMP)/include -I$(NTL)/include $<

.c.o: ecm.h ecm-gmp.h
	$(CC) $(CFLAGS) -I$(GMP)/include -c $<

clean:
	rm -f ecm ecm_with_ntl tune $(OBJS)

dist: $(DIST)
	mkdir ecm-$(VERSION)
	cp $(DIST) $(EXTRADIST) ecm-$(VERSION)
	tar cf ecm-$(VERSION).tar ecm-$(VERSION)
	/bin/rm -fr ecm-$(VERSION)
	gzip --best ecm-$(VERSION).tar

median: median-aux.o auxi.o listz.o toomcook.o median.o
	gcc -lgmp -lm $^ -o median
