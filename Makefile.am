
bin_PROGRAMS = ecm tune

ecm_SOURCES = auxi.c b1_ainc.c bestd.c candi.c ecm.c ecm2.c eval.c getprime.c\
	      listz.c lucas.c main.c pm1.c pp1.c stage2.c toomcook.c trial.c\
	      memory.c mpmod.c mul_lo.c polyeval.c resume.c median.c\
	      smartprp.c schoen_strass.c ks-multiply.c rho.c

tune_SOURCES = mpmod.c tune.c auxi.c mul_lo.c

noinst_HEADERS = ecm.h ecm-gmp.h

EXTRA_DIST = test.pm1 test.pp1 test.ecm ecm.1

man_MANS = ecm.1

AM_LDFLAGS = -lm

check: ecm
	$(srcdir)/test.pp1 ./ecm
	$(srcdir)/test.pm1 ./ecm
	$(srcdir)/test.ecm ./ecm

dist-hook:
	xsltproc -o ecm.1 $(XSLDIR)/manpages/docbook.xsl ecm.xml
