
bin_PROGRAMS = ecm

ecm_SOURCES = auxi.c b1_ainc.c bestd.c candi.c ecm.c ecm2.c eval.c getprime.c\
	      listz.c lucas.c main.c pm1.c pp1.c stage2.c toomcook.c trial.c\
	      memory.c mpmod.c mul_lo.c polyeval.c resume.c median.c\
	      smartprp.c schoen_strass.c ks-multiply.c rho.c

noinst_HEADERS = ecm.h ecm-gmp.h

EXTRA_DIST = test.pm1 test.pp1 test.ecm

AM_LDFLAGS = -lgmp -lm

check: ecm
	$(srcdir)/test.pp1 ./ecm
	$(srcdir)/test.pm1 ./ecm
	$(srcdir)/test.ecm ./ecm
