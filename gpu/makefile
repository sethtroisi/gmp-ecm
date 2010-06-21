stage1 : stage1.o prototype.o getprime.o
	gcc stage1.o prototype.o getprime.o -o stage1 -lgmp -lm
stage1.o : stage1.c prototype.h getprime.h
	gcc -c stage1.c -o stage1.o -lgmp -lm
prototype.o : prototype.c prototype.h
	gcc -c prototype.c -o prototype.o -lgmp -lm
getprime.o : getprime.c getprime.h
	gcc -c getprime.c -o getprime.o -lgmp -lm
