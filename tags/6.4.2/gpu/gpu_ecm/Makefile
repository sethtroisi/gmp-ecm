BIN:=gpu_ecm

# put here the PATH to GMP
GMP_PATH:=/opt/cluster/gmplib
# put here the path to CUDA
CUDA_PATH=/opt/cluster/nvidia/cuda-4.1.28/cuda

# CUDA does not support gcc >= 4.6
# this can be useful if your default gcc is >= 4.6
# if not just comment the following two lines.
#CC_BIN:=/tmp/gcc45
#USE_THIS_CC:=--compiler-bindir $(CC_BIN) 

# Change here the compute capability (can be 2.0 or 1.3)
COMPUTE_CAPABILITY:=2.0

# (can be -m64 or -m32)
mOPTION:=-m64

OBJDIR:=obj

OBJ_CU:=$(addprefix $(OBJDIR)/,$(patsubst %.cu,%.cu.o,$(wildcard *.cu)))
OBJ_C:=$(addprefix $(OBJDIR)/,$(patsubst %.c,%.c.o,$(wildcard *.c)))
FILES_FROM_LIBECM:=random.c batch.c getprime.c
FILES_FROM_ECMECM:=getprime.c eval.c candi.c auxi.c resume.c
OBJ_FROM_LIBECM:=$(addprefix $(OBJDIR)/,$(patsubst %.c,%.l,$(FILES_FROM_LIBECM)))
OBJ_FROM_ECMECM:=$(addprefix $(OBJDIR)/,$(patsubst %.c,%.o,$(FILES_FROM_ECMECM)))
OBJ_FILES:=$(OBJ_C) $(OBJ_CU) $(OBJ_FROM_LIBECM) $(OBJ_FROM_ECMECM) 

LDFLAGS:=$(GMP_PATH)/lib/libgmp.a -lm -lgmp
LDFLAGS+=-L$(CUDA_PATH)/lib64 -L$(CUDA_PATH)/lib -lcudart

CFLAGS:=-O2 -I$(GMP_PATH)/include

NVCCFLAGS:=$(USE_THIS_CC) --compiler-options -fno-strict-aliasing 
NVCCFLAGS+=$(mOPTION) --ptxas-options=-v

ifeq ($(COMPUTE_CAPABILITY), 2.0)
  COMPUTECAPABILITY:=-arch=sm_20
  CFLAGS+=-DCC20
  #NVCCFLAGS+=--maxrregcount 20
endif
ifeq ($(COMPUTE_CAPABILITY), 1.3)
  COMPUTECAPABILITY:=-arch=sm_13
  CFLAGS+=-DCC13
  NVCCFLAGS+=--opencc-options -O2 --maxrregcount 32
endif

#Only for debug
KEEP:=-keep -keep-dir=keep

$(BIN): $(OBJ_FILES)
	gcc -Wall -fPIC $(mOPTION) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.cu.o : %.cu cudakernel.h def.h | $(OBJDIR)
	nvcc $(KEEP) $(CFLAGS) $(COMPUTECAPABILITY) $(NVCCFLAGS) -o $@ -c $<

$(OBJDIR)/%.c.o : %.c def.h utils.h | $(OBJDIR)
	gcc -Wall $(CFLAGS) -o $@ -c $<

$(OBJDIR)/%.o : ../../%.c | $(OBJDIR)
	gcc -Wall $(CFLAGS) -DGPUECM -DVERSION_GPU=\"0.1\" -DOUTSIDE_LIBECM -o $@ -c $<

$(OBJDIR)/%.l : ../../%.c | $(OBJDIR)
	gcc -Wall $(CFLAGS) -DGPUECM -DVERSION_GPU=\"0.1\" -o $@ -c $<

$(OBJDIR): 
	mkdir $(OBJDIR)

clean: 
	rm -f $(OBJDIR)/*.o

distclean: clean
	rm -f $(BIN)
