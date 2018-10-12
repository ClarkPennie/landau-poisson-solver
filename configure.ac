#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([FULL-PACKAGE-NAME], [VERSION], [BUG-REPORT-ADDRESS])
AC_CONFIG_MACRO_DIR([m4])
#AC_CONFIG_SRCDIR([source/collisionRoutines_1.h])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CXX
#AC_PROG_CC
AC_LANG([C++])
ACX_MPI([CC="$MPIXX"], AC_MSG_ERROR([Could not find MPI C++ compiler support]))
CXX="$MPICXX"

# Checks for libraries.
#AX_FFTW3(fftw3)

# simple fftw check

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

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor pow sqrt])

AC_OUTPUT(Makefile source/Makefile)
