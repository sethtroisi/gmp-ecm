
bin_PROGRAMS = ecm

ecm_SOURCES = auxi.c b1_ainc.c bestd.c candi.c ecm.c ecm2.c eval.c getprime.c\
			  listz.c lucas.c main.c median.c memory.c mpmod.c mul_lo.c\
			  pm1.c polyeval.c polyz.c pp1.c resume.c ks-multiply.c\
			  schoen_strass.c smartprp.c stage2.c toomcook.c trial.c rho.c

noinst_HEADERS = ecm.h ecm-gmp.h
EXTRA_DIST = test.pm1 test.pp1 test.ecm

LDFLAGS = -lgmp -lm

test: ecm
	./test.pp1 ./ecm
	./test.pm1 ./ecm
	./test.ecm ./ecm
