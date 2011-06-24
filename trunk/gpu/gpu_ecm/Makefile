BIN=gpu_ecm

GMP=/usr

LDFLAGS=$(GMP)/lib/libgmp.a -lm
CFLAGS=-O2 -I$(GMP)/include
NVCCFLAGS=-m64 --compiler-options -fno-strict-aliasing --opencc-options -O2 --ptxas-options=-v

ifeq ($(keep),1)
NVCCFLAGS+=-keep -keep-dir keep2
endif


ifeq ($(m),1)
CFLAGS+=-DMTEST
endif

ifeq ($(v),0)
CFLAGS+=-DTEST
CFLAGS+=-DSIZE_NUMBER=32
else
	ifeq ($(v),1)
		CFLAGS+=-DSIZE_NUMBER=4
	else
		ifeq ($(v),2)
			CFLAGS+=-DSIZE_NUMBER=8
		else
			ifeq ($(v),3)
				CFLAGS+=-DSIZE_NUMBER=16
			else
				CFLAGS+=-DSIZE_NUMBER=32
			endif
		endif
	endif
endif

ifeq ($(cc),2)
COMPUTECAPABILITY=-arch=sm_20
CFLAGS+=-DCC20
NVCCFLAGS+=--maxrregcount 20
else
COMPUTECAPABILITY=-arch=sm_13
CFLAGS+=-DCC13
NVCCFLAGS+=--maxrregcount 32
endif

ifeq ($(debug),0)
CFLAGS+= 
else
CFLAGS+=-g#-G
endif

$(BIN): obj/utils.cu.o obj/main.cu.o obj/cudautils.cu.o obj/cudaarith.cu.o
	gcc -fPIC -m64 -o $@ $^ $(LDFLAGS) -L/usr/local/cuda/lib64 -lcudart

obj/utils.cu.o : utils.cu utils.h 
	nvcc $(CFLAGS) $(COMPUTECAPABILITY) $(NVCCFLAGS) -o $@ -c $<

obj/main.cu.o : main.cu main.h utils.h cudautils.h
	nvcc $(CFLAGS) $(COMPUTECAPABILITY) $(NVCCFLAGS) -o $@ -c $<

obj/cudautils.cu.o : cudautils.cu main.h cudaarith.h
	nvcc $(CFLAGS) $(COMPUTECAPABILITY) $(NVCCFLAGS) -o $@ -c $<

obj/cudaarith.cu.o : cudaarith.cu main.h
	nvcc $(CFLAGS) $(COMPUTECAPABILITY) $(NVCCFLAGS) -o $@ -c $<

s: s.c
	gcc -I/usr/include -c s.c
	gcc -o stest s.o $(LDFLAGS)

#withdouble: obj/withdouble.cu.o objdir
#obj/withdouble.cu.o : withdouble.cu 
#	nvcc $(CFLAGS) -arch=sm_20 $(NVCCFLAGS) -o $@ -c $<

clean: 
	rm obj/*.cu.o
	
mrproper: clean
	rm $(BIN)

