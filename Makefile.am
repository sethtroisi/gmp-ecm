lib_LIBRARIES = libecm.a

libecm_a_SOURCES = ecm.c ecm2.c pm1.c pp1.c getprime.c listz.c lucas.c\
		   stage2.c toomcook.c mpmod.c mul_lo.c polyeval.c median.c\
		   schoen_strass.c ks-multiply.c rho.c bestd.c auxlib.c\
		   random.c factor.c
if ENABLE_ASM_REDC
  libecm_a_SOURCES += redc.s
endif

if ENABLE_NTT
  libecm_a_SOURCES += sp.c spv.c spm.c mpzspm.c mpzspv.c ntt_gfp.c ecm_ntt.c
endif

if WITH_GWNUM
  libecm_a_SOURCES += Fgw.c
  Fgwtest : Fgw.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -g -DTESTDRIVE -o Fgwtest Fgw.c libecm.a $(LIBS)
endif

bin_PROGRAMS = ecm
noinst_PROGRAMS = tune ecmfactor

LDADD = libecm.a

ecm_SOURCES = auxi.c b1_ainc.c candi.c eval.c random2.c main.c trial.c \
              memory.c resume.c smartprp.c getprime2.c

ecmfactor_SOURCES = ecmfactor.c

tune_CPPFLAGS = -DTUNE $(AM_CPPFLAGS)
tune_LDADD =
tune_SOURCES = mpmod.c tune.c mul_lo.c listz.c auxlib.c ks-multiply.c \
               toomcook.c schoen_strass.c polyeval.c median.c
if ENABLE_NTT
  tune_SOURCES += ecm_ntt.c ntt_gfp.c mpzspv.c mpzspm.c sp.c \
                  spv.c spm.c random.c
endif

include_HEADERS = ecm.h
noinst_HEADERS = ecm-impl.h ecm-gmp.h ecm-ecm.h sp.h longlong.h ecm-params.h

EXTRA_DIST = test.pm1 test.pp1 test.ecm ecm.1 README.lib ecm.xml

ecm-params:	tune
		@echo Optimising parameters for your system, please be patient.
		test -z "ecm-params.h" || rm -f ecm-params.h
		./tune > ecm-params.h

check: ecm
	$(srcdir)/test.pp1 ./ecm
	$(srcdir)/test.pm1 ./ecm
	$(srcdir)/test.ecm ./ecm

man_MANS =
if MAKE_MANPAGE
  man_MANS += ecm.1
endif

ecm.1: ecm.xml
	xsltproc -o ecm.1 $(XSLDIR)/manpages/docbook.xsl ecm.xml
