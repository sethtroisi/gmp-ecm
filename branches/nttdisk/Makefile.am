## Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I m4

# to not install libecm.la, we could write noinst_LTLIBRARIES instead of
# lib_LTLIBRARIES below, however then libecm.a is not installed either
# (see http://www.gnu.org/software/automake/manual/html_node/Libtool-Convenience-Libraries.html)
lib_LTLIBRARIES = libecm.la

EXTRA_PROGRAMS = rho batch

# If we want assembly mulredc code, recurse into the right subdirectory
# and set up variables to include the mulredc library from that subdir
if ENABLE_ASM_REDC
  SUBDIRS = $(ASMPATH)
  MULREDCINCPATH = -I$(srcdir)/$(ASMPATH)
  MULREDCLIBRARY = $(builddir)/$(ASMPATH)/libmulredc.la

  # Add a testing and timing program for the mulredc code
  EXTRA_PROGRAMS += bench_mulredc test_mulredc
  CLEANFILES = bench_mulredc test_mulredc
  bench_mulredc_CPPFLAGS = $(MULREDCINCPATH)
  bench_mulredc_LDADD = $(MULREDCLIBRARY) $(GMPLIB)
  test_mulredc_CPPFLAGS = $(MULREDCINCPATH)
  test_mulredc_LDADD = $(MULREDCLIBRARY) $(GMPLIB)
endif

libecm_la_SOURCES = ecm.c ecm2.c pm1.c pp1.c getprime.c listz.c lucas.c \
		   stage2.c toomcook.c mpmod.c mul_lo.c polyeval.c median.c \
		   schoen_strass.c ks-multiply.c rho.c bestd.c auxlib.c \
		   random.c factor.c sp.c spv.c spm.c mpzspm.c mpzspv.c \
		   ntt_gfp.c ecm_ntt.c pm1fs2.c pm1fs2_param.c mul_fft.c \
		   sets_long.c auxarith.c batch.c
# Link the asm redc code (if we use it) into libecm.la
libecm_la_CPPFLAGS = $(MULREDCINCPATH)
libecm_la_CFLAGS = $(OPENMP_CFLAGS)
libecm_la_LDFLAGS = '-version-info 0:0:0'
libecm_la_LIBADD = $(MULREDCLIBRARY) $(GMPLIB)

bin_PROGRAMS = ecm
noinst_PROGRAMS = tune ecmfactor ecmfactor2
# Most binaries want to link libecm.la, and the ones which don't will
# have their own _LDADD which overrides the default LDADD here
LDADD = libecm.la $(GMPLIB)

ecm_CPPFLAGS = -DOUTSIDE_LIBECM
ecm_CFLAGS = $(OPENMP_CFLAGS)
ecm_SOURCES = auxi.c b1_ainc.c candi.c eval.c random.c main.c trial.c \
              resume.c getprime.c champions.h

tune_SOURCES = mpmod.c tune.c mul_lo.c listz.c auxlib.c ks-multiply.c \
               toomcook.c schoen_strass.c polyeval.c median.c ecm_ntt.c \
	       ntt_gfp.c mpzspv.c mpzspm.c sp.c spv.c spm.c random.c \
	       mul_fft.c auxarith.c
tune_CPPFLAGS = -DTUNE $(MULREDCINCPATH)
tune_LDADD = $(MULREDCLIBRARY) $(GMPLIB)

ecmfactor_CFLAGS = $(OPENMP_CFLAGS)
ecmfactor2_CFLAGS = $(OPENMP_CFLAGS)

rho_SOURCES = rho.c
rho_CPPFLAGS = -DTESTDRIVE
rho_LDADD = -lprimegen $(GMPLIB) $(GSL_LD_FLAGS)

batch_SOURCES = batch.c
batch_LDADD = $(GMPLIB)

if MEMORY_DEBUG
  libecm_la_SOURCES += memory.c
  tune_SOURCES += memory.c
endif

if WITH_GWNUM
  gwdata.ld :
	echo "SECTIONS { .data : { . = ALIGN(0x20); *(_GWDATA) } }" >gwdata.ld
  
  libecm_la_SOURCES += Fgw.c
  # Use ecm_DEPENDENCIES += gwdata.ld instead? Is that possible?
  ecm_DEPENDENCIES = gwdata.ld
  ecm_LDFLAGS = $(AM_LDFLAGS) -Wl,gwdata.ld
  Fgwtest : Fgw.c gwdata.ld
	$(CC) $(CFLAGS) $(CPPFLAGS) -g -DTESTDRIVE -Wl,gwdata.ld -o Fgwtest Fgw.c libecm.a $(LIBS)
endif

include_HEADERS = ecm.h
noinst_HEADERS = ecm-impl.h ecm-gmp.h ecm-ecm.h sp.h longlong.h ecm-params.h

EXTRA_DIST = test.pm1 test.pp1 test.ecm README.lib ecm.xml \
             ecm-params.h.alpha-ev5 ecm-params.h.athlon64 \
             ecm-params.h.default ecm-params.h.alpha-ev6 \
             ecm-params.h.athlon ecm-params.h.powerpc7450 \
             ecm-params.h.pentium3 ecm-params.h.pentium4 \
             ecm-params.h.pentium-m ecm-params.h.powerpc970 \
             ecm-params.h.mips64el ecm-params.h.armv5tel \
             ecm-params.h.sparc64 ecm-params.h.ia64 \
             ecm-params.h.hppa2.0 ecm-params.h.alpha-ev56 \
             ecm-params.h.core2 \
             mul_fft-params.h.athlon64 mul_fft-params.h.pentium3 \
             mul_fft-params.h.default mul_fft-params.h.pentium4

DIST_SUBDIRS = athlon pentium4 x86_64 powerpc64 build.vc9

DISTCLEANFILES = config.m4

ecm-params:	tune$(EXEEXT)
		@echo Optimising parameters for your system, please be patient.
		test -z "ecm-params.h" || rm -f ecm-params.h
		./tune > ecm-params.h

check: ecm$(EXEEXT)
	$(srcdir)/test.pp1 ./ecm$(EXEEXT)
	$(srcdir)/test.pm1 ./ecm$(EXEEXT)
	$(srcdir)/test.ecm ./ecm$(EXEEXT)

longcheck: ecm$(EXEEXT)
	$(srcdir)/test.pp1 "$(VALGRIND) ./ecm$(EXEEXT)"
	$(srcdir)/test.pp1 "$(VALGRIND) ./ecm$(EXEEXT) -no-ntt"
	$(srcdir)/test.pp1 "$(VALGRIND) ./ecm$(EXEEXT) -modmuln"
	$(srcdir)/test.pp1 "$(VALGRIND) ./ecm$(EXEEXT) -redc"
	$(srcdir)/test.pp1 "$(VALGRIND) ./ecm$(EXEEXT) -mpzmod"
	$(srcdir)/test.pm1 "$(VALGRIND) ./ecm$(EXEEXT)"
	$(srcdir)/test.pm1 "$(VALGRIND) ./ecm$(EXEEXT) -no-ntt"
	$(srcdir)/test.pm1 "$(VALGRIND) ./ecm$(EXEEXT) -modmuln"
	$(srcdir)/test.pm1 "$(VALGRIND) ./ecm$(EXEEXT) -redc"
	$(srcdir)/test.pm1 "$(VALGRIND) ./ecm$(EXEEXT) -mpzmod"
	$(srcdir)/test.ecm "$(VALGRIND) ./ecm$(EXEEXT)"
	$(srcdir)/test.ecm "$(VALGRIND) ./ecm$(EXEEXT) -no-ntt"
	$(srcdir)/test.ecm "$(VALGRIND) ./ecm$(EXEEXT) -modmuln"
	$(srcdir)/test.ecm "$(VALGRIND) ./ecm$(EXEEXT) -redc"
	$(srcdir)/test.ecm "$(VALGRIND) ./ecm$(EXEEXT) -mpzmod"

## to be sure ecm.1 is considered as source
## (cf section "Man pages" in the automake manual)
dist_man_MANS = ecm.1

## If ./configure found xsltproc and docbook.xsl, add a rule for building
## the manpage. If they were not found, this rule is omitted and make will
## never try to rebuild the man page, which would result in an error anyway.
if MAKE_MANPAGE
ecm.1: ecm.xml
	xsltproc -o ecm.1 $(XSLDIR)/manpages/docbook.xsl ecm.xml
endif
