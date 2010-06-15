stage1 : stage1.o prototype.o
	gcc stage1.o prototype.o -o stage1 -lgmp -lm
stage1.o : stage1.c prototype.h
	gcc -c stage1.c -o stage1.o -lgmp
prototype.o : prototype.c prototype.h
	gcc -c prototype.c -o prototype.o -lgmp -lm

