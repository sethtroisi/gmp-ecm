## Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = subdir-objects

# to disable valgrind during make distcheck
AM_DISTCHECK_CONFIGURE_FLAGS ?= \
  --disable-valgrind-tests

@COV_FRAG@

# to not install libecm.la, we could write noinst_LTLIBRARIES instead of
# lib_LTLIBRARIES below, however then libecm.a is not installed either
# (see http://www.gnu.org/software/automake/manual/html_node/Libtool-Convenience-Libraries.html)
lib_LTLIBRARIES = libecm.la

EXTRA_PROGRAMS = rho

# If we want assembly mulredc code, recurse into the right subdirectory
# and set up variables to include the mulredc library from that subdir
if ENABLE_ASM_REDC
  SUBDIRS = $(ASMPATH)
  MULREDCINCPATH = -I$(srcdir)/$(ASMPATH)
  MULREDCLIBRARY = $(builddir)/$(ASMPATH)/libmulredc.la

  # Add a tuning and testing program for the mulredc code
  EXTRA_PROGRAMS += bench_mulredc test_mulredc
  CLEANFILES = bench_mulredc test_mulredc
  bench_mulredc_CPPFLAGS = $(MULREDCINCPATH)
  bench_mulredc_LDADD = $(MULREDCLIBRARY) $(GMPLIB)
  test_mulredc_CPPFLAGS = $(MULREDCINCPATH)
  test_mulredc_LDADD = $(MULREDCLIBRARY) $(GMPLIB)
else
  # Add a tuning program for the mulredc code
  EXTRA_PROGRAMS += bench_mulredc
  CLEANFILES = bench_mulredc
  bench_mulredc_LDADD = $(GMPLIB)
endif

libecm_la_SOURCES = ecm.c ecm2.c pm1.c pp1.c getprime_r.c listz.c lucas.c \
		   stage2.c mpmod.c mul_lo.c polyeval.c median.c \
		   schoen_strass.c ks-multiply.c rho.c bestd.c auxlib.c \
		   random.c factor.c sp.c spv.c spm.c mpzspm.c mpzspv.c \
		   ntt_gfp.c ecm_ntt.c pm1fs2.c sets_long.c \
		   auxarith.c batch.c parametrizations.c cudawrapper.c \
		   aprtcle/mpz_aprcl.c addlaws.c torsions.c
# Link the asm redc code (if we use it) into libecm.la
libecm_la_CPPFLAGS = $(MULREDCINCPATH)
libecm_la_CFLAGS = $(OPENMP_CFLAGS) -g
# www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html
# If any interfaces have been added, removed, or changed since the last
# update, increment current, and set revision to 0.
libecm_la_LDFLAGS = $(LIBECM_LDFLAGS) -version-info 1:0:0 -g
libecm_la_LIBADD = $(MULREDCLIBRARY)
if WANT_GPU 
  libecm_la_SOURCES += cudakernel.cu
  libecm_la_LIBADD += $(CUDALIB)
  libecm_la_LDFLAGS += $(CUDALDFLAGS)
  ecm_LDFLAGS = $(CUDARPATH)
endif
libecm_la_LIBADD += $(GMPLIB)

bin_PROGRAMS = ecm
noinst_PROGRAMS = tune ecmfactor bench_mulredc aprcl
# Most binaries want to link libecm.la, and the ones which don't will
# have their own _LDADD which overrides the default LDADD here
LDADD = libecm.la $(GMPLIB)

ecm_CPPFLAGS = -DOUTSIDE_LIBECM
ecm_CFLAGS = $(OPENMP_CFLAGS) -g
ecm_SOURCES = auxi.c b1_ainc.c candi.c eval.c main.c resume.c \
	      addlaws.c torsions.c \
              getprime_r.c champions.h aprtcle/mpz_aprcl.c memusage.c

tune_SOURCES = mpmod.c tune.c mul_lo.c listz.c auxlib.c ks-multiply.c \
               schoen_strass.c polyeval.c median.c ecm_ntt.c \
	       ntt_gfp.c mpzspv.c mpzspm.c sp.c spv.c spm.c auxarith.c
tune_CPPFLAGS = -DTUNE $(MULREDCINCPATH)
tune_LDADD = $(MULREDCLIBRARY) $(GMPLIB)

ecmfactor_CFLAGS = $(OPENMP_CFLAGS) -lpthread

# multiecm_CPPFLAGS = -DOUTSIDE_LIBECM
# multiecm_SOURCES = auxi.c b1_ainc.c candi.c eval.c random.c multiecm.c \
#               resume.c getprime_r.c addlaws.c torsions.c aprtcle/mpz_aprcl.c
# multiecm_CFLAGS = $(OPENMP_CFLAGS)

aprcl_CFLAGS = -Wall
aprcl_SOURCES = aprtcle/mpz_aprcl.c aprtcle/aprcl.c
aprcl_LDADD = $(GMPLIB)

rho_SOURCES = rho.c
rho_CPPFLAGS = -DTESTDRIVE
rho_LDADD = -lprimesieve -lgsl $(GMPLIB)

if WITH_GWNUM
  gwdata.ld :
	echo "SECTIONS { .data : { . = ALIGN(0x20); *(_GWDATA) } }" >gwdata.ld
  
  libecm_la_SOURCES += Fgw.c
  # Use ecm_DEPENDENCIES += gwdata.ld instead? Is that possible?
  ecm_DEPENDENCIES = gwdata.ld
  ecm_LDFLAGS = $(AM_LDFLAGS) -Wl,gwdata.ld
  LIBS += -ldl
  Fgwtest : Fgw.c gwdata.ld
	$(CC) $(CFLAGS) $(CPPFLAGS) -g -DTESTDRIVE -Wl,gwdata.ld -o Fgwtest Fgw.c libecm.a $(LIBS)
endif

include_HEADERS = ecm.h
noinst_HEADERS = basicdefs.h ecm-impl.h ecm-gmp.h ecm-ecm.h sp.h longlong.h \
                 ecm-params.h mpmod.h ecm-gpu.h cudakernel.h addlaws.h \
                 getprime_r.h ecm_int.h \
                 aprtcle/mpz_aprcl.h aprtcle/jacobi_sum.h

EXTRA_DIST = test.pm1 test.pp1 test.ecm README.lib INSTALL-ecm ecm.xml  \
             pentium4/params.h x86_64/core2/params.h x86_64/k8/params.h \
             athlon/params.h x86/params.h ia64/params.h arm/params.h    \
             powerpc64/params.h powerpc32/params.h sparc64/params.h     \
             hppa/params.h mips/params.h x86_64/corei7/params.h         \
             generic/params.h testlong.pp1 testlong.pm1 testlong.ecm \
             getprime_r.h cudakernel_default.cu test.gpuecm README.gpu c155 \
             test_dummy2.save test_dummy.save test_prime95.save M877.save \
	     M997.save test_Z2102.n

DIST_SUBDIRS = athlon pentium4 x86_64 powerpc64 aprtcle build.vc12

DISTCLEANFILES = config.m4

if WANT_GPU
.cu.lo: cudakernel.h
		$(LIBTOOL) --tag=CC --mode=compile $(NVCC) --compile $(NVCCFLAGS) --compiler-options -fPIC $(CPPFLAGS) -o $@ $^ -static
endif

ecm-params:	bench_mulredc$(EXEEXT) tune$(EXEEXT)
		@echo Optimising parameters for your system, please be patient.
		test -z "ecm-params.h" || rm -f ecm-params.h
		./bench_mulredc > ecm-params.h
		./tune >> ecm-params.h

check_PROGRAMS = ecm$(EXEEXT)

dist_check_SCRIPTS = test.pp1 test.pm1 test.ecm
if WANT_GPU
dist_check_SCRIPTS += test.gpuecm
endif

TESTS = $(dist_check_SCRIPTS)
TESTS_ENVIRONMENT = $(VALGRIND)

# see https://www.gnu.org/software/automake/manual/html_node/Scripts_002dbased-Testsuites.html
# this is only needed for "make longcheck" below
export GMPECM_DATADIR = $(srcdir)

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
	$(srcdir)/test.ecm "$(VALGRIND) ./ecm$(EXEEXT) -treefile tree"
	$(srcdir)/testlong.pp1 "$(VALGRIND) ./ecm$(EXEEXT)"
	$(srcdir)/testlong.pm1 "$(VALGRIND) ./ecm$(EXEEXT)"
	$(srcdir)/testlong.ecm "$(VALGRIND) ./ecm$(EXEEXT)"

## to be sure ecm.1 is considered as source
## (cf section "Man pages" in the automake manual)
dist_man_MANS = ecm.1

## If ./configure found xsltproc and docbook.xsl, add a rule for building
## the manpage. If they were not found, this rule is omitted and make will
## never try to rebuild the man page, which would result in an error anyway.
if MAKE_MANPAGE
ecm.1: $(srcdir)/ecm.xml
	xsltproc -o ecm.1 $(XSLDIR)/manpages/docbook.xsl $(srcdir)/ecm.xml
endif
