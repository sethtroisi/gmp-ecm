## Process this file with automake to produce Makefile.in

lib_LIBRARIES = libecm.a

libecm_a_SOURCES = ecm.c ecm2.c pm1.c pp1.c getprime.c listz.c lucas.c \
		   stage2.c toomcook.c mpmod.c mul_lo.c polyeval.c median.c \
		   schoen_strass.c ks-multiply.c rho.c bestd.c auxlib.c \
		   random.c factor.c sp.c spv.c spm.c mpzspm.c mpzspv.c \
		   ntt_gfp.c ecm_ntt.c pm1fs2.c

bin_PROGRAMS = ecm
noinst_PROGRAMS = tune ecmfactor

LDADD = libecm.a

ecm_CPPFLAGS = -DOUTSIDE_LIBECM $(AM_CPPFLAGS)
ecm_SOURCES = auxi.c b1_ainc.c candi.c eval.c random.c main.c trial.c \
              resume.c getprime.c

tune_CPPFLAGS = -DTUNE $(AM_CPPFLAGS)
tune_LDADD =
tune_SOURCES = mpmod.c tune.c mul_lo.c listz.c auxlib.c ks-multiply.c \
               toomcook.c schoen_strass.c polyeval.c median.c ecm_ntt.c \
	       ntt_gfp.c mpzspv.c mpzspm.c sp.c spv.c spm.c random.c

if MEMORY_DEBUG
  libecm_a_SOURCES += memory.c
  tune_SOURCES += memory.c
endif

if ENABLE_ASM_REDC
  nodist_libecm_a_SOURCES = redc.asm mulredc1.asm mulredc2.asm mulredc3.asm \
    mulredc4.asm mulredc5.asm mulredc6.asm mulredc7.asm mulredc8.asm \
    mulredc9.asm mulredc10.asm mulredc11.asm mulredc12.asm mulredc13.asm \
    mulredc14.asm mulredc15.asm mulredc16.asm mulredc17.asm mulredc18.asm \
    mulredc19.asm mulredc20.asm
  nodist_tune_SOURCES = redc.asm mulredc1.asm mulredc2.asm mulredc3.asm \
    mulredc4.asm mulredc5.asm mulredc6.asm mulredc7.asm mulredc8.asm \
    mulredc9.asm mulredc10.asm mulredc11.asm mulredc12.asm mulredc13.asm \
    mulredc14.asm mulredc15.asm mulredc16.asm mulredc17.asm mulredc18.asm \
    mulredc19.asm mulredc20.asm
endif

if WITH_GWNUM
  gwdata.ld :
	echo "SECTIONS { .data : { . = ALIGN(0x20); *(_GWDATA) } }" >gwdata.ld
  
  libecm_a_SOURCES += Fgw.c
  ecm_DEPENDENCIES = $(LDADD) gwdata.ld
  ecm_LDFLAGS = -Wl,gwdata.ld
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
             configfsf.sub configfsf.guess

DIST_SUBDIRS = athlon pentium4 x86_64

ecm-params:	tune
		@echo Optimising parameters for your system, please be patient.
		test -z "ecm-params.h" || rm -f ecm-params.h
		./tune > ecm-params.h

check: ecm
	$(srcdir)/test.pp1 ./ecm
	$(srcdir)/test.pm1 ./ecm
	$(srcdir)/test.ecm ./ecm

.asm.o:
	$(M4) -DOPERATION_$* `test -f $< || echo '$(srcdir)/'`$< >tmp-$*.s
	$(CCAS) $(CFLAGS) tmp-$*.s -o $@
	rm -f tmp-$*.s

## to be sure ecm.1 is considered as source
## (cf section "Man pages" in the automake manual)
dist_man_MANS = ecm.1

man_MANS =
if MAKE_MANPAGE
  man_MANS += ecm.1
endif

ecm.1: ecm.xml
	xsltproc -o ecm.1 $(XSLDIR)/manpages/docbook.xsl ecm.xml

pm1fs2: pm1fs2.c libecm.a
	gcc -Wall -Wextra -DTESTDRIVE -g -O0 -o pm1fs2 pm1fs2.c libecm.a -lm -lgmp
