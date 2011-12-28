

CC=gcc
CFLAGS=-O2 -g -Wall -std=c99




hecm: auxi.o morphismes.o Jacobi.o generation.o ariKS.o hecm.o stage1HECM.o stage2HECM.o
	gcc auxi.o stage1HECM.o morphismes.o Jacobi.o generation.o ariKS.o hecm.o stage2HECM.o ../libecm_la-mpmod.o ../libecm_la-getprime.o ../libecm_la-mul_fft.o ../libecm_la-auxlib.o ../libecm_la-mul_lo.o -o hecm ../.libs/libecm.a -lgmp -lm ../x86_64/.libs/libmulredc.a

generation.o: ../ecm-impl.h generation.h ariKS.h hecm.h
ariKS.o: ../ecm-impl.h ariKS.h generation.h
auxi.o: auxi.h 
morphimes.o: ../ecm-impl.h generation.h ariKS.h morphismes.h Jacobi.h auxi.h 
stage1HECM.o: ../ecm-impl.h morphismes.h ariKS.h generation.h hecm.h
stage2HECM.o: ../ecm-impl.h hecm.h
hecm.o: ../ecm-impl.h morphismes.h generation.h auxi.h hecm.h
Jacobi.o: ../ecm-impl.h generation.h Jacobi.h



clean: 
	rm -f *.o *~ *.bak hecm
