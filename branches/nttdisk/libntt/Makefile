CC32 = gcc -m32
CC64 = gcc -m64

CFLAGS = -O3 -fomit-frame-pointer -I. -Intt
#CFLAGS = -g -I. -Intt

HDR = \
	basicdefs.h \
	gmp-xface.h \
	libntt.h \
	ntt/longlong.h \
	ntt/ntt-impl-scalar.h \
	ntt/ntt-impl-simd.h \
	ntt/ntt-impl-sse2.h \
	ntt/ntt-impl-sse42.h \
	ntt/ntt-impl-avx.h \
	ntt/ntt-impl-avx2.h \
	ntt/ntt-impl-fma.h \
	ntt/ntt-impl.h \
	ntt/sp.h

NTT_SCALAR_SRC = \
	core/ntt02.c \
	core/ntt03.c \
	core/ntt04.c \
	core/ntt05.c \
	core/ntt07.c \
	core/ntt08.c \
	core/ntt09.c \
	core/ntt15.c \
	core/ntt16.c \
	core/ntt35.c \
	core/ntt40.c \
	ntt/mpzspm.c \
	ntt/ntt.c \
	ntt/ntt-scalar.c \
	ntt/sp.c \
	ntt/spm.c \
	ntt/spv.c

NTT_SIMD_SRC = \
	core/ntt02simd.c \
	core/ntt03simd.c \
	core/ntt04simd.c \
	core/ntt05simd.c \
	core/ntt07simd.c \
	core/ntt08simd.c \
	core/ntt09simd.c \
	core/ntt15simd.c \
	core/ntt16simd.c \
	core/ntt35simd.c \
	core/ntt40simd.c \
	ntt/ntt-simd.c

LIB_SRC = \
	libntt.c \
	util.c

OBJ_W32 = $(LIB_SRC:.c=.o_w32) \
	$(NTT_SCALAR_SRC:.c=.o_sp50w32) \
	$(NTT_SIMD_SRC:.c=.o_sp50w32avx) \
	$(NTT_SIMD_SRC:.c=.o_sp50w32fma) \
	$(NTT_SCALAR_SRC:.c=.o_sp30w32) \
	$(NTT_SIMD_SRC:.c=.o_sp30w32sse2) \
	$(NTT_SIMD_SRC:.c=.o_sp30w32sse42) \
	$(NTT_SIMD_SRC:.c=.o_sp30w32avx2) \
	$(NTT_SCALAR_SRC:.c=.o_sp31w32) \
	$(NTT_SIMD_SRC:.c=.o_sp31w32sse2) \
	$(NTT_SIMD_SRC:.c=.o_sp31w32sse42) \
	$(NTT_SIMD_SRC:.c=.o_sp31w32avx2) \
	$(NTT_SCALAR_SRC:.c=.o_sp62w32) \
	$(NTT_SIMD_SRC:.c=.o_sp62w32sse2) \
	$(NTT_SIMD_SRC:.c=.o_sp62w32sse42) \
	$(NTT_SIMD_SRC:.c=.o_sp62w32avx2)

OBJ_W64 = $(LIB_SRC:.c=.o_w64) \
	$(NTT_SCALAR_SRC:.c=.o_sp50w64) \
	$(NTT_SIMD_SRC:.c=.o_sp50w64avx) \
	$(NTT_SIMD_SRC:.c=.o_sp50w64fma) \
	$(NTT_SCALAR_SRC:.c=.o_sp30w64) \
	$(NTT_SIMD_SRC:.c=.o_sp30w64sse2) \
	$(NTT_SIMD_SRC:.c=.o_sp30w64sse42) \
	$(NTT_SIMD_SRC:.c=.o_sp30w64avx2) \
	$(NTT_SCALAR_SRC:.c=.o_sp31w64) \
	$(NTT_SIMD_SRC:.c=.o_sp31w64sse2) \
	$(NTT_SIMD_SRC:.c=.o_sp31w64sse42) \
	$(NTT_SIMD_SRC:.c=.o_sp31w64avx2) \
	$(NTT_SCALAR_SRC:.c=.o_sp62w64) \
	$(NTT_SIMD_SRC:.c=.o_sp62w64sse2) \
	$(NTT_SIMD_SRC:.c=.o_sp62w64sse42) \
	$(NTT_SIMD_SRC:.c=.o_sp62w64avx2)

w64: $(OBJ_W64)
	$(CC64) $(CFLAGS) $(OBJ_W64) test.c -o test64 -lgmp -lm

w32: $(OBJ_W32)
	$(CC32) $(CFLAGS) $(OBJ_W32) test.c -o test32 -lgmp -lm

clean:
	rm -rf $(OBJ_W32) $(OBJ_W64) test32.exe test64.exe ./test

%.o_w32: %.c $(HDR)
	$(CC32) $(CFLAGS) -c -o $@ $<

%.o_sp30w32: %.c $(HDR)
	$(CC32) $(CFLAGS) -DHAVE_SSE2 -DHAVE_SSE42 -DHAVE_AVX2 -DSP_NUMB_BITS=30 -c -o $@ $<

%.o_sp30w32sse2: %.c $(HDR)
	$(CC32) $(CFLAGS) -march=core2 -DHAVE_SSE2 -DSP_NUMB_BITS=30 -c -o $@ $<

%.o_sp30w32sse42: %.c $(HDR)
	$(CC32) $(CFLAGS) -march=nehalem -DHAVE_SSE42 -DSP_NUMB_BITS=30 -c -o $@ $<

%.o_sp30w32avx2: %.c $(HDR)
	$(CC32) $(CFLAGS) -march=haswell -DHAVE_AVX2 -DSP_NUMB_BITS=30 -c -o $@ $<

%.o_sp31w32: %.c $(HDR)
	$(CC32) $(CFLAGS) -DHAVE_SSE2 -DHAVE_SSE42 -DHAVE_AVX2 -DSP_NUMB_BITS=31 -c -o $@ $<

%.o_sp31w32sse2: %.c $(HDR)
	$(CC32) $(CFLAGS) -march=core2 -DHAVE_SSE2 -DSP_NUMB_BITS=31 -c -o $@ $<

%.o_sp31w32sse42: %.c $(HDR)
	$(CC32) $(CFLAGS) -march=nehalem -DHAVE_SSE42 -DSP_NUMB_BITS=31 -c -o $@ $<

%.o_sp31w32avx2: %.c $(HDR)
	$(CC32) $(CFLAGS) -march=haswell -DHAVE_AVX2 -DSP_NUMB_BITS=31 -c -o $@ $<

%.o_sp62w32: %.c $(HDR)
	$(CC32) $(CFLAGS) -DHAVE_SSE2 -DHAVE_SSE42 -DHAVE_AVX2 -DSP_NUMB_BITS=62 -c -o $@ $<

%.o_sp62w32sse2: %.c $(HDR)
	$(CC32) $(CFLAGS) -march=core2 -DHAVE_SSE2 -DSP_NUMB_BITS=62 -c -o $@ $<

%.o_sp62w32sse42: %.c $(HDR)
	$(CC32) $(CFLAGS) -march=nehalem -DHAVE_SSE42 -DSP_NUMB_BITS=62 -c -o $@ $<

%.o_sp62w32avx2: %.c $(HDR)
	$(CC32) $(CFLAGS) -march=haswell -DHAVE_AVX2 -DSP_NUMB_BITS=62 -c -o $@ $<

%.o_sp50w32: %.c $(HDR)
	$(CC32) $(CFLAGS) -DHAVE_AVX -DHAVE_FMA -DSP_NUMB_BITS=50 -c -o $@ $<

%.o_sp50w32avx: %.c $(HDR)
	$(CC32) $(CFLAGS) -march=ivybridge -DHAVE_AVX -DSP_NUMB_BITS=50 -c -o $@ $<

%.o_sp50w32fma: %.c $(HDR)
	$(CC32) $(CFLAGS) -march=haswell -DHAVE_FMA -DSP_NUMB_BITS=50 -c -o $@ $<


%.o_w64: %.c $(HDR)
	$(CC64) $(CFLAGS) -c -o $@ $<

%.o_sp30w64: %.c $(HDR)
	$(CC64) $(CFLAGS) -DHAVE_SSE2 -DHAVE_SSE42 -DHAVE_AVX2 -DSP_NUMB_BITS=30 -c -o $@ $<

%.o_sp30w64sse2: %.c $(HDR)
	$(CC64) $(CFLAGS) -march=core2 -DHAVE_SSE2 -DSP_NUMB_BITS=30 -c -o $@ $<

%.o_sp30w64sse42: %.c $(HDR)
	$(CC64) $(CFLAGS) -march=nehalem -DHAVE_SSE42 -DSP_NUMB_BITS=30 -c -o $@ $<

%.o_sp30w64avx2: %.c $(HDR)
	$(CC64) $(CFLAGS) -march=haswell -DHAVE_AVX2 -DSP_NUMB_BITS=30 -c -o $@ $<

%.o_sp31w64: %.c $(HDR)
	$(CC64) $(CFLAGS) -DHAVE_SSE2 -DHAVE_SSE42 -DHAVE_AVX2 -DSP_NUMB_BITS=31 -c -o $@ $<

%.o_sp31w64sse2: %.c $(HDR)
	$(CC64) $(CFLAGS) -march=core2 -DHAVE_SSE2 -DSP_NUMB_BITS=31 -c -o $@ $<

%.o_sp31w64sse42: %.c $(HDR)
	$(CC64) $(CFLAGS) -march=nehalem -DHAVE_SSE42 -DSP_NUMB_BITS=31 -c -o $@ $<

%.o_sp31w64avx2: %.c $(HDR)
	$(CC64) $(CFLAGS) -march=haswell -DHAVE_AVX2 -DSP_NUMB_BITS=31 -c -o $@ $<

%.o_sp62w64: %.c $(HDR)
	$(CC64) $(CFLAGS) -DHAVE_SSE2 -DHAVE_SSE42 -DHAVE_AVX2 -DSP_NUMB_BITS=62 -c -o $@ $<

%.o_sp62w64sse2: %.c $(HDR)
	$(CC64) $(CFLAGS) -march=core2 -DHAVE_SSE2 -DSP_NUMB_BITS=62 -c -o $@ $<

%.o_sp62w64sse42: %.c $(HDR)
	$(CC64) $(CFLAGS) -march=nehalem -DHAVE_SSE42 -DSP_NUMB_BITS=62 -c -o $@ $<

%.o_sp62w64avx2: %.c $(HDR)
	$(CC64) $(CFLAGS) -march=haswell -DHAVE_AVX2 -DSP_NUMB_BITS=62 -c -o $@ $<

%.o_sp50w64: %.c $(HDR)
	$(CC64) $(CFLAGS) -DHAVE_AVX -DHAVE_FMA -DSP_NUMB_BITS=50 -c -o $@ $<

%.o_sp50w64avx: %.c $(HDR)
	$(CC64) $(CFLAGS) -march=ivybridge -DHAVE_AVX -DSP_NUMB_BITS=50 -c -o $@ $<

%.o_sp50w64fma: %.c $(HDR)
	$(CC64) $(CFLAGS) -march=haswell -DHAVE_FMA -DSP_NUMB_BITS=50 -c -o $@ $<

