
lib_LIBRARIES = libecm.a

libecm_a_SOURCES = ecm.c ecm2.c pm1.c pp1.c getprime.c listz.c lucas.c\
		   stage2.c toomcook.c mpmod.c mul_lo.c polyeval.c median.c\
		   schoen_strass.c ks-multiply.c rho.c bestd.c auxlib.c\
		   random.c factor.c

bin_PROGRAMS = ecm
noinst_PROGRAMS = tune tune2 ecmfactor

ecm_SOURCES = auxi.c b1_ainc.c candi.c eval.c random2.c\
	      main.c trial.c memory.c resume.c smartprp.c getprime2.c
ecmfactor_SOURCES = ecmfactor.c

tune_SOURCES = mpmod.c tune.c mul_lo.c
tune2_SOURCES = tune2.c

LDADD = libecm.a

include_HEADERS = ecm.h
noinst_HEADERS = ecm-impl.h ecm-gmp.h ecm-ecm.h

EXTRA_DIST = test.pm1 test.pp1 test.ecm ecm.1 README.lib ecm.xml

man_MANS = ecm.1

check: ecm
	$(srcdir)/test.pp1 ./ecm
	$(srcdir)/test.pm1 ./ecm
	$(srcdir)/test.ecm ./ecm

ecm.1: ecm.xml
	xsltproc -o ecm.1 $(XSLDIR)/manpages/docbook.xsl ecm.xml
