BIN=gpu_ecm

# put here the PATH to GMP
GMP_PATH=/usr
# put here the path to CUDA
CUDA_PATH=/usr/local/cuda

# CUDA does not support gcc >= 4.6
# this can be useful if your default gcc is >= 4.6
# if not just comment the following two lines.
CC_BIN=/tmp/gcc45
USE_THIS_CC=--compiler-bindir $(CC_BIN) 

# Change here the compute capability (can be 2.0 or 1.3)
COMPUTE_CAPABILITY=1.3

# (can be -m64 or -m32)
mOPTION=-m64

OBJDIR=obj
#CU_FILES=cudakernel.cu main.cu cudautils.cu cudaarith.cu
#C_FILES=utils.c 
OBJ_FILES=$(addprefix $(OBJDIR)/,utils.c.o main.c.o cudakernel.cu.o)

LDFLAGS=$(GMP_PATH)/lib/libgmp.a -lm 
LDFLAGS+=-L$(CUDA_PATH)/lib64 -lcudart

CFLAGS=-O2 -I$(GMP_PATH)/include -I$(CUDA_PATH)/include

NVCCFLAGS=$(USE_THIS_CC) --compiler-options -fno-strict-aliasing 
NVCCFLAGS+=$(mOPTION) --ptxas-options=-v

ifeq ($(COMPUTE_CAPABILITY), 2.0)
  COMPUTECAPABILITY=-arch=sm_20
  CFLAGS+=-DCC20
  NVCCFLAGS+=--maxrregcount 20
endif
ifeq ($(COMPUTE_CAPABILITY), 1.3)
  COMPUTECAPABILITY=-arch=sm_13
  CFLAGS+=-DCC13
  NVCCFLAGS+=--opencc-options -O2 --maxrregcount 32
endif


$(BIN): $(OBJ_FILES) | $(OBJDIR)
	gcc -Wall -fPIC $(mOPTION) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.cu.o : %.cu
	nvcc $(CFLAGS) $(COMPUTECAPABILITY) $(NVCCFLAGS) -o $@ -c $<

$(OBJDIR)/%.c.o : %.c
	gcc -Wall $(CFLAGS) -o $@ -c $<

$(OBJDIR): 
	mkdir $(OBJDIR)

clean: 
	rm -f $(OBJ_FILES)

distclean: clean
	rm -f $(BIN)
