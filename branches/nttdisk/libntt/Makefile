CC32 = mingw32-gcc
CC64 = gcc -I../../../mpir-2.7.0/lib/x64/Release

CFLAGS = -O3 -march=native -fomit-frame-pointer -I. -Intt

HDR = \
	basicdefs.h \
	gmp-xface.h \
	libntt.h \
	ntt/longlong.h \
	ntt/ntt-impl-sse2.h \
	ntt/ntt-impl.h \
	ntt/sp.h

CORE_SRC = \
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
	core/ntt40.c

DRIVER_SRC = \
	ntt/mpzspm.c \
	ntt/ntt.c \
	ntt/sp.c \
	ntt/spm.c \
	ntt/spv.c \
	ntt/test.c

NTT_SRC = $(CORE_SRC) $(DRIVER_SRC)

LIBNTT_SRC = \
	util.c

OBJ_W32 = $(LIBNTT_SRC:.c=.o_w32)
OBJ_W64 = $(LIBNTT_SRC:.c=.o_w64)

OBJ_SP30_W32 = $(NTT_SRC:.c=.o_sp30w32)
OBJ_SP30_W32_SSE2 = $(NTT_SRC:.c=.o_sp30w32sse2)
OBJ_SP31_W32 = $(NTT_SRC:.c=.o_sp31w32)
OBJ_SP31_W32_SSE2 = $(NTT_SRC:.c=.o_sp31w32sse2)
OBJ_SP62_W32 = $(NTT_SRC:.c=.o_sp62w32)
OBJ_SP62_W32_SSE2 = $(NTT_SRC:.c=.o_sp62w32sse2)

OBJ_SP30_W64 = $(NTT_SRC:.c=.o_sp30w64)
OBJ_SP30_W64_SSE2 = $(NTT_SRC:.c=.o_sp30w64sse2)
OBJ_SP62_W64 = $(NTT_SRC:.c=.o_sp62w64)
OBJ_SP62_W64_SSE2 = $(NTT_SRC:.c=.o_sp62w64sse2)

OBJ_ALL = \
	$(OBJ_W32) \
	$(OBJ_W64) \
	$(OBJ_SP30_W32) \
	$(OBJ_SP30_W32_SSE2) \
	$(OBJ_SP31_W32) \
	$(OBJ_SP31_W32_SSE2) \
	$(OBJ_SP62_W32) \
	$(OBJ_SP62_W32_SSE2) \
	$(OBJ_SP30_W64) \
	$(OBJ_SP30_W64_SSE2) \
	$(OBJ_SP62_W64) \
	$(OBJ_SP62_W64_SSE2)

w32: $(OBJ_W32) $(OBJ_SP30_W32) $(OBJ_SP30_W32_SSE2) \
		$(OBJ_SP31_W32) $(OBJ_SP31_W32_SSE2) \
		$(OBJ_SP62_W32) $(OBJ_SP62_W32_SSE2)
	$(CC32) $(CFLAGS) $(OBJ_W32) $(OBJ_SP30_W32_SSE2) -o a32.exe -lgmp

w64: $(OBJ_W64) $(OBJ_SP30_W64) $(OBJ_SP30_W32_SSE2) \
		$(OBJ_SP62_W64) $(OBJ_SP62_W64_SSE2)
	$(CC64) $(CFLAGS) $(OBJ_W64) $(OBJ_SP62_W64_SSE2) -o a64.exe \
		../../../mpir-2.7.0/lib/x64/Release/mpir.lib \
		../../../mpir-2.7.0/lib/x64/Release/chkstk.obj  \
		../../../mpir-2.7.0/lib/x64/Release/runtmchk.lib

clean:
	rm -rf $(OBJ_ALL)

%.o_w32: %.c $(HDR)
	$(CC32) $(CFLAGS) -c -o $@ $<

%.o_sp30w32: %.c $(HDR)
	$(CC32) $(CFLAGS) -UHAVE_SSE2 -DSP_NUMB_BITS=30 -c -o $@ $<

%.o_sp30w32sse2: %.c $(HDR)
	$(CC32) $(CFLAGS) -DHAVE_SSE2 -DSP_NUMB_BITS=30 -c -o $@ $<

%.o_sp31w32: %.c $(HDR)
	$(CC32) $(CFLAGS) -UHAVE_SSE2 -DSP_NUMB_BITS=31 -c -o $@ $<

%.o_sp31w32sse2: %.c $(HDR)
	$(CC32) $(CFLAGS) -DHAVE_SSE2 -DSP_NUMB_BITS=31 -c -o $@ $<

%.o_sp62w32: %.c $(HDR)
	$(CC32) $(CFLAGS) -UHAVE_SSE2 -DSP_NUMB_BITS=62 -c -o $@ $<

%.o_sp62w32sse2: %.c $(HDR)
	$(CC32) $(CFLAGS) -DHAVE_SSE2 -DSP_NUMB_BITS=62 -c -o $@ $<


%.o_w64: %.c $(HDR)
	$(CC64) $(CFLAGS) -c -o $@ $<

%.o_sp30w64: %.c $(HDR)
	$(CC64) $(CFLAGS) -UHAVE_SSE2 -DSP_NUMB_BITS=30 -c -o $@ $<

%.o_sp30w64sse2: %.c $(HDR)
	$(CC64) $(CFLAGS) -DHAVE_SSE2 -DSP_NUMB_BITS=30 -c -o $@ $<

%.o_sp62w64: %.c $(HDR)
	$(CC64) $(CFLAGS) -UHAVE_SSE2 -DSP_NUMB_BITS=62 -c -o $@ $<

%.o_sp62w64sse2: %.c $(HDR)
	$(CC64) $(CFLAGS) -DHAVE_SSE2 -DSP_NUMB_BITS=62 -c -o $@ $<

