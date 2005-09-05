## Process this file with automake to produce Makefile.in

lib_LIBRARIES = libecm.a

libecm_a_SOURCES = ecm.c ecm2.c pm1.c pp1.c getprime.c listz.c lucas.c \
		   stage2.c toomcook.c mpmod.c mul_lo.c polyeval.c median.c \
		   schoen_strass.c ks-multiply.c rho.c bestd.c auxlib.c \
		   random.c factor.c sp.c spv.c spm.c mpzspm.c mpzspv.c \
		   ntt_gfp.c ecm_ntt.c

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
  libecm_a_SOURCES += redc.asm
  tune_SOURCES += redc.asm
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
noinst_HEADERS = ecm-impl.h ecm-gmp.h ecm-ecm.h sp.h longlong.h ecm-params.h

EXTRA_DIST = test.pm1 test.pp1 test.ecm ecm.1 README.lib ecm.xml \
             ecm-params.h.alpha-ev5 ecm-params.h.athlon64 \
             ecm-params.h.default ecm-params.h.alpha-ev6 \
             ecm-params.h.athlon ecm-params.h.powerpc7450 \
             configfsf.sub configfsf.guess


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
	$(CCAS) $(COMPILE_FLAGS) tmp-$*.s -o $@
	rm -f tmp-$*.s


man_MANS =
if MAKE_MANPAGE
  man_MANS += ecm.1
endif

ecm.1: ecm.xml
	xsltproc -o ecm.1 $(XSLDIR)/manpages/docbook.xsl ecm.xml
