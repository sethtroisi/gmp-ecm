BIN=gpu_ecm

GMP=/usr

LDFLAGS=$(GMP)/lib/libgmp.a -lm
CFLAGS=-O2 -I$(GMP)/include
NVCCFLAGS=-m64 --compiler-options -fno-strict-aliasing --opencc-options -O2 --ptxas-options=-v

COMPUTECAPABILITY=-arch=sm_13
CFLAGS+=-DCC13


$(BIN): obj/main.cu.o obj/cudautils.cu.o obj/cudaarith.cu.o
	gcc -fPIC -m64 -o $@ $^ $(LDFLAGS) -L/usr/local/cuda/lib64 -lcudart

obj/main.cu.o : main.cu cudautils.cu utils.h objdir 
	nvcc $(CFLAGS) $(COMPUTECAPABILITY) $(NVCCFLAGS) -o $@ -c $<

obj/cudautils.cu.o : cudautils.cu cudaarith.cu utils.h objdir
	nvcc $(CFLAGS) $(COMPUTECAPABILITY) $(NVCCFLAGS) -o $@ -c $<

obj/cudaarith.cu.o : cudaarith.cu utils.h objdir
	nvcc $(CFLAGS) $(COMPUTECAPABILITY) $(NVCCFLAGS) -o $@ -c $<

objdir:
	mkdir -p obj

clean: 
	rm obj/*.cu.o
	
mrproper: clean
	rm $(BIN)
