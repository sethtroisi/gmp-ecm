## Process this file with automake to produce Makefile.in

lib_LTLIBRARIES = libecm.la

libecm_la_SOURCES = ecm.c ecm2.c pm1.c pp1.c getprime.c listz.c lucas.c \
		   stage2.c toomcook.c mpmod.c mul_lo.c polyeval.c median.c \
		   schoen_strass.c ks-multiply.c rho.c bestd.c auxlib.c \
		   random.c factor.c sp.c spv.c spm.c mpzspm.c mpzspv.c \
		   ntt_gfp.c ecm_ntt.c pm1fs2.c mul_fft.c sets_long.c \
		   auxarith.c

libecm_la_LDFLAGS = '-version-info 0:0:0'

bin_PROGRAMS = ecm
noinst_PROGRAMS = tune ecmfactor

LDADD = libecm.la

ecm_CPPFLAGS = -DOUTSIDE_LIBECM $(AM_CPPFLAGS)
ecm_SOURCES = auxi.c b1_ainc.c candi.c eval.c random.c main.c trial.c \
              resume.c getprime.c
# We add the list of speed critical object files to the ecm linker flags 
# which isn't very clean, but it seems to be the only way of getting them
# in front of the full list of input object files. Automake does not allow
# modifying the ecm_OBJECTS variable directly.
ecm_LDFLAGS = 
ecm_LDFLAGS += $(ECM_FASTOBJ)

tune_CPPFLAGS = -DTUNE $(AM_CPPFLAGS)
# tune_LDADD =
tune_SOURCES = mpmod.c tune.c mul_lo.c listz.c auxlib.c ks-multiply.c \
               toomcook.c schoen_strass.c polyeval.c median.c ecm_ntt.c \
	       ntt_gfp.c mpzspv.c mpzspm.c sp.c spv.c spm.c random.c \
	       mul_fft.c auxarith.c

if MEMORY_DEBUG
  libecm_la_SOURCES += memory.c
  tune_SOURCES += memory.c
endif

if ENABLE_ASM_REDC
# The code in redc.S needs preprocessing by cpp, since we need different
# addressing in PIC and non-PIC mode. Hence the upper-case .S, which
# signifies a file that is to be run through cpp before assembly
  nodist_libecm_la_SOURCES = redc.S mulredc1.s mulredc2.s mulredc3.s \
    mulredc4.s mulredc5.s mulredc6.s mulredc7.s mulredc8.s \
    mulredc9.s mulredc10.s mulredc11.s mulredc12.s mulredc13.s \
    mulredc14.s mulredc15.s mulredc16.s mulredc17.s mulredc18.s \
    mulredc19.s mulredc20.s
# .s and .S files don't get cleaned by default?
  CLEANFILES = redc.S mulredc1.s mulredc2.s mulredc3.s \
    mulredc4.s mulredc5.s mulredc6.s mulredc7.s mulredc8.s \
    mulredc9.s mulredc10.s mulredc11.s mulredc12.s mulredc13.s \
    mulredc14.s mulredc15.s mulredc16.s mulredc17.s mulredc18.s \
    mulredc19.s mulredc20.s
  tune_LDADD = redc.lo mulredc1.o mulredc2.o mulredc3.o \
    mulredc4.o mulredc5.o mulredc6.o mulredc7.o mulredc8.o \
    mulredc9.o mulredc10.o mulredc11.o mulredc12.o mulredc13.o \
    mulredc14.o mulredc15.o mulredc16.o mulredc17.o mulredc18.o \
    mulredc19.o mulredc20.o
endif

if WITH_GWNUM
  gwdata.ld :
	echo "SECTIONS { .data : { . = ALIGN(0x20); *(_GWDATA) } }" >gwdata.ld
  
  libecm_la_SOURCES += Fgw.c
  # Use ecm_DEPENDENCIES += gwdata.ld instead? Is that possible?
  ecm_DEPENDENCIES = $(LDADD) gwdata.ld
  ecm_LDFLAGS += -Wl,gwdata.ld
  Fgwtest : Fgw.c gwdata.ld
	$(CC) $(CFLAGS) $(CPPFLAGS) -g -DTESTDRIVE -Wl,gwdata.ld -o Fgwtest Fgw.c libecm.a $(LIBS)
endif

ecmfactor_SOURCES = ecmfactor.c

include_HEADERS = ecm.h
noinst_HEADERS = ecm-impl.h ecm-gmp.h ecm-ecm.h sp.h longlong.h ecm-params.h \
                 asmredc.h

EXTRA_DIST = test.pm1 test.pp1 test.ecm README.lib ecm.xml \
             ecm-params.h.alpha-ev5 ecm-params.h.athlon64 \
             ecm-params.h.default ecm-params.h.alpha-ev6 \
             ecm-params.h.athlon ecm-params.h.powerpc7450 \
             ecm-params.h.pentium3 ecm-params.h.pentium4 \
             ecm-params.h.pentium-m configfsf.sub configfsf.guess \
             mul_fft-params.h.athlon64 mul_fft-params.h.pentium3 \
             mul_fft-params.h.default mul_fft-params.h.pentium4 \
             mul_fft-params.h


DIST_SUBDIRS = athlon pentium4 x86_64 build.vc9

DISTCLEANFILES = config.m4

ecm-params:	tune$(EXEEXT)
		@echo Optimising parameters for your system, please be patient.
		test -z "ecm-params.h" || rm -f ecm-params.h
		./tune > ecm-params.h

check: ecm$(EXEEXT)
	$(srcdir)/test.pp1 ./ecm$(EXEEXT)
	$(srcdir)/test.pm1 ./ecm$(EXEEXT)
	$(srcdir)/test.ecm ./ecm$(EXEEXT)

.asm.s:
	$(M4) -DOPERATION_$* `test -f $< || echo '$(srcdir)/'`$< >$*.s
.asm.S:
	$(M4) -DOPERATION_$* `test -f $< || echo '$(srcdir)/'`$< >$*.S

## to be sure ecm.1 is considered as source
## (cf section "Man pages" in the automake manual)
dist_man_MANS = ecm.1

man_MANS =
if MAKE_MANPAGE
  man_MANS += ecm.1
endif

ecm.1: ecm.xml
	xsltproc -o ecm.1 $(XSLDIR)/manpages/docbook.xsl ecm.xml
