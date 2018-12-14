AC_PREREQ([2.69])
AC_INIT([landau], [0.2.0], [BUG-REPORT-ADDRESS])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([/source/collisionRoutines_1.h])
AC_CONFIG_HEADERS([config.h])
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE
AX_SUMMARIZE_ENV

# Check for MPI toolchain

AC_PROG_CXX
AC_LANG([C])
ACX_MPI([CC="$MPICC"], AC_MSG_ERROR([Could not find MPI C compiler support]))
AC_LANG([C++])
ACX_MPI([CXX="$MPICXX"], AC_MSG_ERROR([Could not find MPI C++ compiler support]))

# Check compiler supports the C++11 standard (required for GRVY)
AX_CXX_COMPILE_STDCXX_11([noext], [mandatory])

# OpenMP check
AX_OPENMP([CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"], AC_MSG_ERROR([Could not detect OpenMP linkage]))

# Math library check(s) - the first preference is for MKL. If MKL not
# found, then check for OpenBLAS

# MKL check
AC_LANG_PUSH([C])
AX_COMPILER_VENDOR
ACX_MKL([],AC_MSG_NOTICE([MKL not found. Will next check for OpenBLAS]))
AC_LANG_POP([C])

if test "x${MKLROOT}" != "x"; then
  AC_MSG_NOTICE([MKL detected])
else
  MKL_LIBS=""
fi

# OpenBLAS check
if test "x${MKLROOT}" == "x"; then
AC_LANG_PUSH([C])
LDFLAGS="-L$OPENBLAS_LIB"
#AX_COMPILER_VENDOR
AX_BLAS([CXXFLAGS="$CXXFLAGS -DHAVE_OPENBLAS"],AC_MSG_ERROR([Could not detect OpenBLAS linkage]))
AC_LANG_POP([C])
fi

# GRVY Check
AX_PATH_GRVY([0.32],[yes])
CPPFLAGS="${GRVY_CFLAGS} ${CPPFLAGS}"
LDFLAGS="${GRVY_LIBS} ${LDFLAGS}"

# A simple fftw check
AC_ARG_WITH(fftw, 
  [AS_HELP_STRING([--with-fftw[=DIR]],[root directory of FFTW installation (default = FFTW_DIR)])],
  [with_fftw=$withval
if test "${with_fftw}" != yes; then
    FFTW_PREFIX=$withval
fi
],[
with_fftw=$withval
if test "x${FFTW_DIR}" != "x"; then
   FFTW_PREFIX=${FFTW_DIR}
else
   if test "x${TACC_FFTW3_DIR}" != "x"; then
      FFTW_PREFIX=${TACC_FFTW3_DIR}
   fi
fi
])

if test -d "${FFTW_PREFIX}/lib" ; then
   FFTW_LIBS="-L${FFTW_PREFIX}/lib -lfftw3_threads -lfftw3"
fi

if test -d "${FFTW_PREFIX}/include" ; then
   FFTW_CFLAGS="-I${FFTW_PREFIX}/include"
fi


ac_FFTW_save_CPPFLAGS="$CPPFLAGS"
ac_FFTW_save_LDFLAGS="$LDFLAGS"
ac_FFTW_save_LIBS="$LIBS"

CPPFLAGS="${FFTW_CFLAGS} ${CPPFLAGS}"
LDFLAGS="${FFTW_LIBS} ${LDFLAGS}"
AC_LANG_PUSH([C++])
AC_CHECK_HEADER([fftw3.h],[found_header=yes],[found_header=no])
AC_CHECK_LIB([fftw3],fftw_init_threads,[found_library=yes],[found_library=no])
AC_LANG_POP([C++])

succeeded=no
if test "$found_header" = yes; then
    if test "$found_library" = yes; then
        succeeded=yes
    fi
fi

if test "$succeeded" = no; then
  AC_MSG_ERROR([FFTW not found.  Try either --with-fftw or setting FFTW_DIR.])
fi

AC_SUBST(FFTW_CFLAGS)
AC_SUBST(FFTW_LIBS)
AC_SUBST(FFTW_PREFIX)

# Checks for header files.
AC_CHECK_HEADERS([malloc.h stdlib.h])

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor pow sqrt])

AC_OUTPUT(Makefile source/Makefile)
AX_SUMMARIZE_CONFIG
