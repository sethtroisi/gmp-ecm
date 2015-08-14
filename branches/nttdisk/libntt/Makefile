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

NTT_SRC = \
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
	ntt/sp.c \
	ntt/spm.c \
	ntt/spv.c \
	ntt/test.c

LIB_SRC = \
	util.c

OBJ_W32 = $(LIB_SRC:.c=.o_w32) \
	$(NTT_SRC:.c=.o_sp30w32) \
	$(NTT_SRC:.c=.o_sp31w32) \
	$(NTT_SRC:.c=.o_sp62w32)

OBJ_W64 = $(LIB_SRC:.c=.o_w64) \
	$(NTT_SRC:.c=.o_sp30w64) \
	$(NTT_SRC:.c=.o_sp62w64)

w64: $(OBJ_W64)
	$(CC64) $(CFLAGS) $(OBJ_W64) test.c -o a64.exe \
		../../../mpir-2.7.0/lib/x64/Release/mpir.lib \
		../../../mpir-2.7.0/lib/x64/Release/chkstk.obj  \
		../../../mpir-2.7.0/lib/x64/Release/runtmchk.lib

w32: $(OBJ_W32)
	$(CC32) $(CFLAGS) $(OBJ_W32) test.c -o a32.exe -lgmp

clean:
	rm -rf $(OBJ_W32) $(OBJ_W64) a32.exe a64.exe

%.o_w32: %.c $(HDR)
	$(CC32) $(CFLAGS) -c -o $@ $<

%.o_sp30w32: %.c $(HDR)
	$(CC32) $(CFLAGS) -DHAVE_SSE2 -DSP_NUMB_BITS=30 -c -o $@ $<

%.o_sp31w32: %.c $(HDR)
	$(CC32) $(CFLAGS) -DHAVE_SSE2 -DSP_NUMB_BITS=31 -c -o $@ $<

%.o_sp62w32: %.c $(HDR)
	$(CC32) $(CFLAGS) -DHAVE_SSE2 -DSP_NUMB_BITS=62 -c -o $@ $<


%.o_w64: %.c $(HDR)
	$(CC64) $(CFLAGS) -c -o $@ $<

%.o_sp30w64: %.c $(HDR)
	$(CC64) $(CFLAGS) -DHAVE_SSE2 -DSP_NUMB_BITS=30 -c -o $@ $<

%.o_sp62w64: %.c $(HDR)
	$(CC64) $(CFLAGS) -DHAVE_SSE2 -DSP_NUMB_BITS=62 -c -o $@ $<

