m4_define([ECM_VERSION_AC], [7.0.5-dev])

AC_PREREQ([2.57])
AC_INIT([ecm], ECM_VERSION_AC, [ecm-discuss@inria.fr])

AC_CONFIG_MACRO_DIR([m4])
GMP_INIT([config.m4])
AM_INIT_AUTOMAKE
AC_CANONICAL_HOST

AC_CONFIG_HEADERS([config.h ecm.h])
AC_DEFINE([ECM_VERSION], ["ECM_VERSION_AC"], [The version of GMP-ECM])

dnl Copied from MPFR 2.4.2:
unset GMP_CFLAGS GMP_CC user_redefine_cc
dnl Check if user request his CC and CFLAGS
if test -n "$CFLAGS" || test -n "$CC" ; then
  user_redefine_cc=yes
fi

dnl the following is required to compile auxi.c according to autoconf 2.61
AC_PROG_EGREP
AC_PROG_SED

AC_ARG_WITH([gmp],         [  --with-gmp=DIR          GMP install directory ],  [with_gmp_include=$withval/include with_gmp_lib=$withval/lib])
AC_ARG_WITH([gmp_include], [  --with-gmp-include=DIR  GMP include directory ],  [with_gmp_include=$withval])
AC_ARG_WITH([gmp_lib],     [  --with-gmp-lib=DIR      GMP lib directory ],      [with_gmp_lib=$withval])
AC_ARG_WITH([gwnum],       [  --with-gwnum=DIR        GWNUM source directory ], [with_gwnum=$withval])

if test x"$with_gmp_include" != "x" 
then
  if ! test -d "$with_gmp_include"
  then
    AC_MSG_ERROR([Specified GMP include directory "$with_gmp_include" does not exist])
  fi
  CPPFLAGS="-I$with_gmp_include"
fi

if test x"$with_gmp_lib" != "x" 
then
  if ! test -d "$with_gmp_lib"
  then
    AC_MSG_ERROR([Specified GMP library directory "$with_gmp_lib" does not exist])
  fi
fi

############################
# Parse --enable-* options #
############################

dnl Assertions are enabled by default for beta/rc releases. The last parameter
dnl of AC_ARG_ENABLE() sets the default value (change also default=...).
AC_ARG_ENABLE([assert],
[AS_HELP_STRING([--enable-assert], [enable ASSERT checking [[default=yes]]])],[],[enable_assert=yes])
if test "x$enable_assert" = xyes; then
  AC_DEFINE([WANT_ASSERT],1,[Define to 1 if you want assertions enabled])
  GMP_DEFINE([WANT_ASSERT], 1)
fi

AC_ARG_ENABLE([gmp-cflags],
[AS_HELP_STRING([--enable-gmp-cflags], [enable importing CFLAGS from gmp.h [[default=yes]]])],[],[enable_gmp_cflags=yes])

AC_ARG_ENABLE([openmp],
[AS_HELP_STRING([--enable-openmp], [enable OpenMP multi-threading [[default=no]]])])

AC_ARG_ENABLE([sse2],
[AS_HELP_STRING([--enable-sse2], [use SSE2 instructions in NTT code (default=yes for 32-bit x86 systems, if supported)])])

AC_ARG_ENABLE([aprcl],
[AS_HELP_STRING([--enable-aprcl], [use APRCL to prove factors prime [[default=yes]]])])

AC_ARG_ENABLE([asm-redc],
[AS_HELP_STRING([--enable-asm-redc], [use an asm redc (default=yes on x86_64 and powerpc64/elf, no on others)])])

AC_ARG_ENABLE([mulredc-svoboda],
[AS_HELP_STRING([--enable-mulredc-svoboda], [enable Svoboda mulredc [[default=no]]])])
if test "x$enable_mulredc_svoboda" = xyes; then
  AC_DEFINE([MULREDC_SVOBODA],1,[Define to 1 if you want Svoboda mulredc])
  GMP_DEFINE([MULREDC_SVOBODA], 1)
fi

AC_ARG_ENABLE([valgrind-client],
[AS_HELP_STRING([--enable-valgrind-client], [enable Valgrind client check requests [[default=no]]])],[],[])
if test "x$enable_valgrind" = xyes; then
  AC_DEFINE([USE_VALGRIND],1,[Define to 1 if you want Valgrind client check requests enabled])
fi


dnl Use C language for test programs
AC_LANG([C])

dnl Copied from MPFR 2.4.1 and modified
dnl We need to guess the C preprocessor instead of using AC_PROG_CPP,
dnl since AC_PROG_CPP implies AC_PROG_CC, which chooses a compiler
dnl (before we have the chance to get it from gmp.h) and does some
dnl checking related to this compiler (such as dependency tracking
dnl options); if the compiler changes due to __GMP_CC in gmp.h, one
dnl would have incorrect settings.
dnl FIXME: Move this in aclocal ?

if test "x$user_redefine_cc" = x && test "x$enable_gmp_cflags" = xyes && test "x$cross_compiling" != xyes; then
if test "x$GMP_CC$GMP_CFLAGS" = x; then
   AC_MSG_CHECKING([for CC and CFLAGS in gmp.h])
   GMP_CC=__GMP_CC
   GMP_CFLAGS=__GMP_CFLAGS
   for cpp in /lib/cpp gcc cc c99
   do
     case $cpp in
       *cpp*) cpp="$cpp -P" ;;
       *) cpp="$cpp -E -P" ;;
     esac
     echo "Trying to run $cpp" >&AS_MESSAGE_LOG_FD
     AC_LANG_CONFTEST([AC_LANG_SOURCE([foo])])
     if $cpp $CPPFLAGS conftest.$ac_ext > /dev/null 2> /dev/null ; then
       # Get CC and CFLAGS
       AC_LANG_CONFTEST([AC_LANG_SOURCE([#include "gmp.h"
MPFR_OPTION_CC __GMP_CC
MPFR_OPTION_CFLAGS __GMP_CFLAGS])])
       echo "Trying to parse gmp.h with: $cpp $CPPFLAGS conftest.$ac_ext" >&AS_MESSAGE_LOG_FD
       if $cpp $CPPFLAGS conftest.$ac_ext 2> /dev/null > conftest.out; then
         GMP_CC="`$EGREP MPFR_OPTION_CC conftest.out | $SED -e 's/MPFR_OPTION_CC //g' | $SED -e 's/"//g'`"
         GMP_CFLAGS="`$EGREP MPFR_OPTION_CFLAGS conftest.out | $SED -e 's/MPFR_OPTION_CFLAGS //g'| $SED -e 's/"//g'`"
         echo "Success, GMP_CC=$GMP_CC, GMP_CFLAGS=$GMP_CFLAGS" >&AS_MESSAGE_LOG_FD
         break
       else
         echo "Could not parse gmp.h with $cpp" >&AS_MESSAGE_LOG_FD
       fi
     else
       echo "Could not run $cpp" >&AS_MESSAGE_LOG_FD
     fi
   done
   rm -f conftest*
   if test "x$GMP_CC" = "x__GMP_CC" || test "x$GMP_CFLAGS" = "x__GMP_CFLAGS" ; then
      AC_MSG_RESULT([no])
      GMP_CFLAGS=
      GMP_CC=
   else
      AC_MSG_RESULT([yes CC=$GMP_CC CFLAGS=$GMP_CFLAGS])
   fi
fi

dnl But these variables may be invalid, so we must check them first.
dnl Note: we do not use AC_RUN_IFELSE, as it implies AC_PROG_CC.
if test "x$GMP_CC$GMP_CFLAGS" != x; then
   AC_MSG_CHECKING([whether CC=$GMP_CC and CFLAGS=$GMP_CFLAGS works])
   AC_LANG_CONFTEST([AC_LANG_SOURCE([int main (void) { return 0; }])])
   if $GMP_CC $GMP_CFLAGS -o conftest conftest.$ac_ext 2> /dev/null ; then
     AC_MSG_RESULT([yes])
     CFLAGS=$GMP_CFLAGS
     CC=$GMP_CC
   else
     AC_MSG_RESULT([no, reverting to default])
   fi
   rm -f conftest*
fi

fi


dnl Checks for programs.
AC_PROG_CC_C99
AM_PROG_AS
AM_PROG_CC_C_O


# Now that we have decided on CC and CFLAGS, init libtool
# Don't make a shared library by default. Enable building a shared library
# by specifying "--enable-shared" on the ./configure command line
LT_PREREQ([2.2.6])
LT_INIT([disable-shared,win32-dll])


dnl Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_INT64_T
AC_TYPE_UINT64_T
AC_TYPE_LONG_LONG_INT

dnl Check if compiler supports "const," if not define it to empty string
AC_C_CONST

dnl Check if compiler supports "inline," if not define it to 
dnl __inline__ or __inline or the empty string
AC_C_INLINE

dnl Check if both time.h and sys/time.h can be included
AC_HEADER_TIME

dnl Define size_t to something useable if standard headers don't
AC_TYPE_SIZE_T

dnl Use "zu" for the size_t format string, and "zd" for ssize_t
AC_DEFINE([PRIuSIZE], "zu", [Define to the format string for a size_t])
AC_DEFINE([PRIdSIZE], "zd", [Define to the format string for an ssize_t])

dnl If OpenMP is enabled, check which command line parameter (if any)
dnl if required to make the compiler enable OpenMP
if test "x$enable_openmp" = xyes; then
  AC_OPENMP
fi

# Determine if assembly code is ELF
AC_MSG_CHECKING([if assembly code is ELF])
AC_EGREP_CPP(yes,
[#if __ELF__
  yes
#endif
],
[elf=yes],
[elf=no])
AC_MSG_RESULT([$elf])

########################
# Enable asm redc code #
########################
# If --(en|dis)able-asm-redc not specified, choose default value
# Test if asm redc code is available for this cpu.
# Point ASMPATH to the correct subdirectory.

# asm_redc enabled by default for x86_64 and 64 bit PowerPC
if test "x$enable_asm_redc" = x; then
  case $host in
    x86_64*-*-* | powerpc-apple-darwin*) enable_asm_redc=yes;;
    # powerpc64 assembly is ELF-specific
    powerpc64-*-linux*) enable_asm_redc=$elf;;
    *) enable_asm_redc=no;;
  esac
fi

if test "x$enable_asm_redc" = xyes; then
  case $host in
    pentium4-*-* | pentium3-*-* | viac7-*-* | i786-*-*) ASMPATH=pentium4;;
    x86_64*-*-*)
# In case GMP has been compiled with a 32-bit ABI...
# Use AC_COMPILE_IFELSE instead of AC_PREPROC_IFELSE, otherwise
# GMP's CFLAGS doesn't seem to be taken into account.
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([
#if defined(__i386__)
#error
#endif])], [], [AC_MSG_NOTICE([32-bit ABI (i386), disabling asm-redc])
                enable_asm_redc=no])
                          ASMPATH=x86_64;;
# warning: with powerpc-apple-darwin* we can have ABI=32
# see bug #10646 on the bug tracker, where config.guess says
# powerpc-apple-darwin8.11.0 (this a 64-bit machine, but most applications
# are compiled in 32 bits). It works with --disable-asm-redc.
    powerpc-apple-darwin*)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([
#if defined(__ppc__)
#error
#endif])], [], [AC_MSG_NOTICE([32-bit PowerPC, disabling asm-redc])
                enable_asm_redc=no])
                          ASMPATH=powerpc64;;
    powerpc64-*-linux*)
			  ECM_INCLUDE([<"$srcdir"/powerpc64/powerpc-defs.m4>])
			  ASMPATH=powerpc64;;
    i[[56]]86-*-* | k[[78]]*-*-* | athlon*-*-* | pentiumpro-*-* | \
    pentium2-*-* | viac3*-*-* | i686-apple-darwin*)
    		   	  ASMPATH=athlon;;
    *) AC_MSG_ERROR([[asm redc not available on this machine $host]]);;
  esac
fi

if test "x$enable_asm_redc" = xyes; then
# do the necessary definitions and includes
  AC_DEFINE([USE_ASM_REDC],1,[Define to 1 to use asm redc])
  test "x$CCAS" != x || CCAS="$CC -c"
  AC_SUBST([CCAS])
  GMP_PROG_M4
  GMP_ASM_UNDERSCORE
  GMP_ASM_TEXT
  GMP_ASM_GLOBL
  GMP_ASM_TYPE 

  case $host in
    *-*-mingw32|*cygwin*) GMP_DEFINE([WINDOWS64_ABI], 1)
                 AC_DEFINE([WINDOWS64_ABI], 1,[Define to 1 if x86_64 mulredc*() functions should be called with Windows ABI]);;
    *) ;;
  esac

  case $host in
    pentium3-*-*)
      echo "WARNING: Your processor is recognized as Pentium3."
      echo "         The asm code uses SSE2, and therefore it might"
      echo "         fail if your proc is indeed a P3, and not a"
      echo "         Pentium M. If you have compilation problems,"
      echo "         consider using --disable-asm-redc." ;;
    *)
  esac
fi
AM_CONDITIONAL([ENABLE_ASM_REDC], [test "x$enable_asm_redc" = xyes])


############################
# Enable SSE2 instructions #
############################
# Test if we should use SSE2 instructions and if the cpu supports them
if test "x$enable_sse2" = "x"; then
dnl Default: if we build for Pentium 4, enable SSE2 code for the NTT
dnl Some SSE2 enabled cpus are identified as i686, we enable SSE2 
dnl for them by default and let the tests below check if it works
  case $host in
    pentium4-*-* | viac7-*-* | i686-*-* | i786-*-*) 
      enable_sse2=yes
    ;;
  esac
fi

if test "x$enable_sse2" = xyes; then
  # See if we need -msse2 to enable SSE2 instructions
  AC_MSG_CHECKING([for SSE2 support with inline assembly])
  ECM_RUN_IFELSE([ECM_C_INLINESSE2_PROG], dnl
   [ dnl SSE2 works, nothing to be done
     AC_MSG_RESULT([yes])], dnl
   [ dnl SSE2 does not work, try again with -msse2
    OLDCFLAGS="$CFLAGS"
    CFLAGS="$CFLAGS -msse2"
    ECM_RUN_IFELSE([ECM_C_INLINESSE2_PROG], dnl
      [ dnl works now, keep CFLAGS like this
        AC_MSG_RESULT([yes, with -msse2])], dnl
      [ dnl still didn't work, switch off SSE2
       CFLAGS="$OLDCFLAGS"
       enable_sse2=no
       AC_MSG_RESULT([not supported, SSE2 disabled])
      ])])
fi
if test "x$enable_sse2" = xyes; then
  AC_DEFINE([HAVE_SSE2],1,[Define to 1 to enable SSE2 instructions in NTT code])
fi

#####################
# Enable aprcl code #
#####################
# If --(en|dis)able-aprcl not specified, choose default value

# aprcl enabled by default for all systems
if test "x$enable_aprcl" = x; then
  enable_aprcl=yes
fi

if test "x$enable_aprcl" = xyes; then
  AC_DEFINE([HAVE_APRCL],1,[Define to 1 to enable APRCL primality proving of found factors])
fi

########################
# Add GWNUM if desired #
########################
dnl If user wants GWNUM, check if the file exists (either as .a or .lib)
if test "x$with_gwnum" != "x"; then
  if test "x$enable_openmp" = xyes; then
      AC_MSG_ERROR([Woltman's GWNUM currently cannot be used together with OpenMP])
  fi
  AC_CHECK_FILE([$with_gwnum/gwnum.a], [
    AC_DEFINE([HAVE_GWNUM], 1, [Define to 1 if gwnum.a or gwnum.lib exist])
    GWLIB="$with_gwnum/gwnum.a -lpthread"
    CPPFLAGS="$CPPFLAGS -I$with_gwnum"
  ],[
    AC_CHECK_FILE([$with_gwnum/gwnum.lib], [
      AC_DEFINE([HAVE_GWNUM], 1, [Define to 1 if gwnum.a or gwnum.lib exist])
      GWLIB="$with_gwnum/gwnum.lib -lpthread"
      CPPFLAGS="$CPPFLAGS -I$with_gwnum"
    ],[
      with_gwnum=
      AC_MSG_ERROR([Woltman's GWNUM library not found])
    ])
  ])
fi
AM_CONDITIONAL([WITH_GWNUM], [test "x$with_gwnum" != "x"])


dnl Checks for header files.
AC_FUNC_ALLOCA
AC_HEADER_STDC
AC_CHECK_HEADERS([math.h limits.h malloc.h strings.h sys/time.h unistd.h io.h signal.h fcntl.h])
AC_CHECK_HEADERS([windows.h psapi.h], [], [],
   [[#ifdef HAVE_WINDOWS_H
     # include <windows.h>
     #endif
     ]])
AC_CHECK_HEADERS([ctype.h sys/types.h sys/resource.h aio.h])

dnl Checks for library functions that are not in GMP
AC_FUNC_STRTOD

dnl Check functions in the math library
AC_CHECK_LIB(m,pow,,AC_MSG_ERROR(required function missing))
AC_CHECK_LIB(m,floor,,AC_MSG_ERROR(required function missing))
AC_CHECK_LIB(m,sqrt,,AC_MSG_ERROR(required function missing))
AC_CHECK_LIB(m,fmod,,AC_MSG_ERROR(required function missing))
AC_CHECK_LIB(m,cos)
AC_CHECK_LIB(rt,aio_read)
AC_CHECK_LIB(psapi,[GetProcessMemoryInfo])

AC_CHECK_FUNCS([isascii memset strchr strlen strncasecmp strstr], [], [AC_MSG_ERROR([required function missing])])
AC_CHECK_FUNCS([access unlink], [], [AC_MSG_ERROR([required function missing])])
AC_CHECK_FUNCS([isspace isdigit isxdigit], [], [AC_MSG_ERROR([required function missing])])
AC_CHECK_FUNCS([time ctime], [], [AC_MSG_ERROR([required function missing])])
AC_CHECK_FUNCS([gethostname gettimeofday getrusage memmove signal fcntl fileno setvbuf fallocate aio_read aio_init])

dnl Test for some Windows-specific functions that are available under MinGW
dnl FIXME: which win32 library contains these functions?
dnl AC_CHECK_FUNCS([GetCurrentProcess GetProcessTimes])
AC_CHECK_FUNCS([_fseeki64 _ftelli64])
AC_CHECK_FUNCS([malloc_usable_size])


dnl If we use GCC and user has not specified his own CFLAGS, 
dnl add some warning flags, avoiding duplication
if test "x$GCC" = xyes && test "x$user_redefine_cc" != xyes; then
  case $CFLAGS in
    "-pedantic "* | *" -pedantic "* | *" -pedantic") ;;
    *) CFLAGS="-pedantic $CFLAGS"
  esac
  case $CFLAGS in
    "-Wundef "* | *" -Wundef "* | *" -Wundef") ;;
    *) CFLAGS="-Wundef $CFLAGS"
  esac
  case $CFLAGS in
    "-Wall "* | *" -Wall "* | *" -Wall") ;;
    *) CFLAGS="-Wall $CFLAGS"
  esac
  case $CFLAGS in
    "-W "* | *" -W "* | *" -W") ;;
    *) CFLAGS="-W $CFLAGS"
  esac
#  CFLAGS="-Wall -W -Wundef -pedantic $CFLAGS"
fi


dnl Find GMP and do some sanity checks
dnl Tests concerning the include directories.
if test -d "$with_gmp_include"; then
  dnl AC_CHECK_HEADERS and AC_PREPROC_IFELSE uses CPPFLAGS but not CFLAGS
  CPPFLAGS="-I$with_gmp_include $CPPFLAGS"
else
  with_gmp_include=
fi
AC_CHECK_HEADERS([gmp.h], [], AC_MSG_ERROR([required header file missing]))

dnl This needs to find only the header file so we can do it here, before
dnl we start looking for libgmp.a
AC_MSG_CHECKING([for recent GMP])
AC_PREPROC_IFELSE([AC_LANG_SOURCE([[
#include <gmp.h>
#if (__GNU_MP_VERSION <= 4)
#error
#IRIXdoesnotexitaterrordirective
#endif
]])],
  [AC_MSG_RESULT(yes)],
  [AC_MSG_RESULT(no)
   AC_MSG_ERROR([GMP 5.0.0 or newer is required])]
)

dnl Figure out if GMP is MPIR
AC_MSG_CHECKING([if GMP is MPIR])
AC_PREPROC_IFELSE([AC_LANG_SOURCE([[
#include <gmp.h>
#ifndef __MPIR_VERSION
#error
#IRIXdoesnotexitaterrordirective
#endif
]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_MPIR],1,[Define to 1 if GMP is MPIR])],
  [AC_MSG_RESULT(no)])

dnl Now choose how to link the GMP library. If we can, we'd prefer to
dnl link it statically since that makes for faster function calls. To
dnl link it statically, we mustn't build a dynamic GMP-ECM library and
dnl we need to find the libgmp.a file. At the moment, we only look for
dnl it at the path specified by the user (i.e. --with-gmp) but make no
dnl attempt to find it in the default system lib directories.

dnl If GMP is linked statically, we pass its path/filename to the Makefile
dnl via GMPLIB, otherwise -lgmp is passed via GMPLIB.
dnl The search path to the dynamic GMP library is added to LDFLAGS, 
dnl if GMP is not specified by full pathname.

GMPLDFLAGS=""
if test -d "$with_gmp_lib"; then
  GMPLDFLAGS="-L$with_gmp_lib"
fi
GMPLIB="-lgmp"
dnl On OS X we should not give libgmp.a
dnl (see https://trac.sagemath.org/ticket/20385#comment:103)
case $host_os in
  darwin*) ;;
  *) if test "x$enable_shared" != xyes; then
       if test -r "$with_gmp_lib/libgmp.a"; then
         GMPLIB="$with_gmp_lib/libgmp.a"
         dnl Don't need -L flag since we give full path to libgmp.a
         GMPLDFLAGS=""
       fi
     fi
esac
AC_SUBST([GMPLIB])
LDFLAGS="$LDFLAGS $GMPLDFLAGS"


dnl Test linking against GMP. This tests, for example, that the compiler
dnl and GMP agree on the ABI (32 or 64 bit). AC_CHECK_LIB() does not work,
dnl as it requires linking the library with -lgmp, whereas we may want
dnl to specify it by full pathname as we do in the Makefile

AC_MSG_CHECKING([whether we can link against GMP])
LIBS_BACKUP="$LIBS"
LIBS="$LIBS $GMPLIB"
AC_LINK_IFELSE([AC_LANG_PROGRAM([[
#ifdef HAVE_GMP_H
#include <gmp.h>
#endif]],[[
  mpz_t t;
  mpz_init(t);
  mpz_clear(t);
  return 0;
]])], [AC_MSG_RESULT(yes)],
   [AC_MSG_RESULT(no)
    AC_MSG_ERROR([Could not link against GMP library.])]
)

dnl Check for corresponding 'gmp.h' and libgmp.a
dnl This requires running a program linked against GMP, 
dnl so is done after the link-only test.
AC_MSG_CHECKING([if gmp.h version and libgmp version are the same])
AC_RUN_IFELSE([AC_LANG_PROGRAM([[
#include <stdio.h>
#include <string.h>
#include <gmp.h>
]],[[
char buffer[100];
if (__GNU_MP_VERSION == 4 && __GNU_MP_VERSION_MINOR <= 2 && 
    __GNU_MP_VERSION_PATCHLEVEL == 0)
  sprintf (buffer, "%d.%d", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR);
else
  sprintf (buffer, "%d.%d.%d", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR,
   __GNU_MP_VERSION_PATCHLEVEL);
printf ("(%s/%s) ", buffer, gmp_version);
return strcmp (buffer, gmp_version);
]])], [AC_MSG_RESULT([yes])],
   [AC_MSG_RESULT([no])
    AC_MSG_ERROR(['gmp.h' and 'libgmp' have different versions, you have to reinstall GMP properly, or use the --with-gmp parameter to tell configure the path to the GMP library and header you want to use])],
   [AC_MSG_RESULT([cross-compiling: cannot test])]
)

AC_CHECK_FUNCS([__gmpn_add_nc __gmpn_mod_34lsub1 __gmpn_redc_1 __gmpn_redc_2])
AC_CHECK_FUNCS([__gmpn_mullo_n __gmpn_redc_n __gmpn_preinv_mod_1 __gmpn_mod_1s_4p_cps __gmpn_mod_1s_4p])
AC_CHECK_FUNCS([__gmpn_mul_fft __gmpn_fft_next_size __gmpn_fft_best_k])
AC_CHECK_FUNCS([__gmpn_mulmod_bnm1 __gmpn_mulmod_bnm1_next_size])

LIBS="$LIBS_BACKUP"


dnl Copied from MPFR 3.1.1 and modified.
dnl Configs for Windows DLLs.
dnl libtool requires "-no-undefined" for win32 dll
dnl
dnl "-Wl,output-def" is used to get a .def file for use by MS lib to make
dnl a .lib import library, described in the manual.
dnl
dnl -version-info is set to 0.0.0 in Makefile.am, so we just use 0 in the
dnl .def file filename as well
dnl
dnl Incidentally, libtool does generate an import library libecm.dll.a,
dnl but it's "ar" format and cannot be used by the MS linker.  There
dnl doesn't seem to be any GNU tool for generating or converting to .lib.
case $host in
  *-*-cygwin* | *-*-mingw* | *-*-pw32* | *-*-os2*)
   AC_MSG_CHECKING(for DLL/static GMP)
   if test "$enable_shared" = yes; then
     LIBECM_LDFLAGS="$LIBECM_LDFLAGS -no-undefined -Wl,--output-def,.libs/libecm-0.def"
     AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include "gmp.h"
#if !__GMP_LIBGMP_DLL
# error "Dead man"
error
#endif
     ]], [[]])],[AC_MSG_RESULT(DLL)],[
  AC_MSG_RESULT(static)
  AC_MSG_ERROR([gmp.h isn't a DLL: use --enable-static --disable-shared]) ])
   else
     AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include "gmp.h"
#if __GMP_LIBGMP_DLL
# error "Dead man"
error
#endif
     ]], [[]])],[AC_MSG_RESULT(static)],[
  AC_MSG_RESULT(DLL)
  AC_MSG_ERROR([gmp.h is a DLL: use --disable-static --enable-shared]) ])
  fi

  AC_MSG_CHECKING([whether we can link without -lpsapi])
  # AC_SEARCH_LIBS won't work here, so do it by hand...
  AC_DEFUN([ECM_PSAPI_PROGRAM],
    [AC_LANG_PROGRAM([[
	#ifdef _WIN32
	#include <windows.h>
	#include <psapi.h>
	int testfunc ()
	{
	    PROCESS_MEMORY_COUNTERS info;
	    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
	    return (int) (info.PeakWorkingSetSize >> 10);
	}
	#endif
	]])])
  AC_LINK_IFELSE([ECM_PSAPI_PROGRAM],
    [AC_MSG_RESULT(yes)],
    [AC_MSG_RESULT(no)
    try_link_psapi=yes
   ])
  if test x"$try_link_psapi" = "xyes"; then
    AC_MSG_CHECKING([whether we can link with -lpsapi])
    LIBS="-lpsapi $LIBS"
    AC_LINK_IFELSE([ECM_PSAPI_PROGRAM],
      [AC_MSG_RESULT(yes)],
      [AC_MSG_RESULT(no)
       AC_MSG_ERROR([memusage.c using unresolved GetProcessMemoryInfo])
     ])
  fi
  ;;
esac
AC_SUBST([LIBECM_LDFLAGS])


dnl Check if the compiler understands some __attribute__ directives
AC_MSG_CHECKING([whether compiler knows __attribute__((hot))])
dnl The AC_LANG_WERROR directive causes configure to consider a test 
dnl compilation unsuccessful if it produced any output on stderr.
dnl We use it since unknown __attribute__ only cause a warning, not an 
dnl error. Unfortunately there seems to be no way to switch it off again, 
dnl so this test is at the end of the configure script
AC_LANG_WERROR
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([void foo() __attribute__ ((hot));
void foo() {return;}], [foo()])],
[AC_DEFINE([ATTRIBUTE_HOT],[__attribute__ ((hot))], [How to specify hot-spot attribute, if available])
 AC_MSG_RESULT([yes])
],
[AC_DEFINE([ATTRIBUTE_HOT],[ ], [How to specify hot-spot attribute, if available])
 AC_MSG_RESULT([no])
])

dnl Check for xsltproc
AC_CHECK_PROG([XSLTPROC],[xsltproc],[xsltproc])
if test "x$XSLTPROC" != x; then
  AC_MSG_CHECKING([for docbook.xsl])
  if test "x$XSLDIR" = x; then
    if test -d "/usr/local/share/sgml/docbook/xsl-stylesheets"; then
      XSLDIR="/usr/local/share/sgml/docbook/xsl-stylesheets"
    elif test -d "/usr/share/sgml/docbook/xsl-stylesheets"; then
      XSLDIR="/usr/share/sgml/docbook/xsl-stylesheets"
    elif test -d "/usr/local/share/docbook/"; then
      XSLDIR="/usr/local/share/docbook/"
    elif test -d "/usr/share/docbook/"; then
      XSLDIR="/usr/share/docbook/"
    fi
  fi
  if test -r "$XSLDIR/manpages/docbook.xsl"; then
    AC_MSG_RESULT([yes])
    make_manpage="yes"
  else
    AC_MSG_RESULT([no])
  fi
fi  
AM_CONDITIONAL([MAKE_MANPAGE], [test "x$make_manpage" = xyes])

dnl Check for valgrind. GMP-ECM uses exit code 1 to signal error,
dnl so we make valgrind use that code to signal error, too
dnl AC_CHECK_PROG([VALGRIND], [valgrind], [valgrind -q --error-exitcode=1])
gl_VALGRIND_TESTS

dnl Check if a GPU version is asked, for which GPU and where CUDA is install.
dnl See acinclude.m4 for details.
CU_CHECK_CUDA

LIBS="$LIBS $GWLIB"

GMP_FINISH
AC_CONFIG_FILES([Makefile athlon/Makefile pentium4/Makefile x86_64/Makefile powerpc64/Makefile aprtcle/Makefile build.vc12/Makefile build.vc12/assembler/Makefile build.vc12/ecm/Makefile build.vc12/ecm_gpu/Makefile build.vc12/libecm/Makefile build.vc12/libecm_gpu/Makefile build.vc12/tune/Makefile build.vc12/bench_mulredc/Makefile])

AC_SUBST([XSLDIR])
AC_SUBST([ASMPATH])

COV_FRAG='
# These options are used for coverage tests
PARAMS00 ?= 0
PARAMS11 ?= 0
PARAMS22 ?= 0
PARAMS33 ?= 0
ifeq (${PARAMS00}, 1)
  DEFS += -DPARAMS00
endif
ifeq (${PARAMS11}, 1)
  DEFS += -DPARAMS11
endif
ifeq (${PARAMS22}, 1)
  DEFS += -DPARAMS22
endif
ifeq (${PARAMS33}, 1)
  DEFS += -DPARAMS33
endif
'
AC_SUBST([COV_FRAG])
AM_SUBST_NOTMAKE([COV_FRAG])

AC_OUTPUT

AC_MSG_NOTICE([Configuration:])
AC_MSG_NOTICE([Build for host type $host])

AC_MSG_NOTICE([CC=$CC, CFLAGS=$CFLAGS])
AC_MSG_NOTICE([Linking GMP with $GMPLIB])

if test "x$enable_asm_redc" = xyes; then
  AC_MSG_NOTICE([Using asm redc code from directory $ASMPATH])
else
  AC_MSG_NOTICE([Not using asm redc code])
fi

if test "x$enable_sse2" = xyes; then
  AC_MSG_NOTICE([Using SSE2 instructions in NTT code])
else
  AC_MSG_NOTICE([Not using SSE2 instructions in NTT code])
fi

if test "x$enable_aprcl" = xyes; then
  AC_MSG_NOTICE([Using APRCL to prove factors prime/composite])
else
  AC_MSG_NOTICE([Not using APRCL])
fi

if test "x$with_gwnum" != "x"; then
  AC_MSG_NOTICE([Linking with George Woltman's GWNUM])
fi

if test "x$enable_assert" = xyes; then
  AC_MSG_NOTICE([Assertions enabled])
else
  AC_MSG_NOTICE([Assertions disabled])
fi

if test "x$enable_shellcmd" = xyes; then
  AC_MSG_NOTICE([Shell command execution enabled])
else
  AC_MSG_NOTICE([Shell command execution disabled])
fi

if test "x$enable_openmp" = xyes; then
  AC_MSG_NOTICE([OpenMP enabled])
else
  AC_MSG_NOTICE([OpenMP disabled])
fi

if test "x$enable_memory_debug" = xyes; then
  AC_MSG_NOTICE([Memory debugging enabled])
else
  AC_MSG_NOTICE([Memory debugging disabled])
fi
