# SYNOPSIS
#
#   Summarizes configuration settings.
#
#   AX_SUMMARIZE_CONFIG([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Outputs a summary of relevant configuration settings.
#
# LAST MODIFICATION
#
#   2018-10-12
#

AC_DEFUN([AX_SUMMARIZE_CONFIG],
[

echo
echo '----------------------------------- SUMMARY -----------------------------------'
echo
echo Package version............... : $PACKAGE-$VERSION
echo
echo C++ compiler.................. : $CXX
echo C++ compiler flags............ : $CXXFLAGS
echo Install dir................... : $prefix 
echo FFTW includes................. : $FFTW_CFLAGS
echo FFTW lib...................... : $FFTW_LIBS
echo MKLROOT....................... : $MKLROOT
echo MKL_LIBS...................... : $MKL_LIBS
echo BLAS_LIBS..................... : $BLAS_LIBS
echo Build user.................... : $USER
echo Build host.................... : $BUILD_HOST
echo Configure date................ : $BUILD_DATE
echo Build architecture............ : $BUILD_ARCH
echo Git revision number........... : $BUILD_VERSION
echo
echo '-------------------------------------------------------------------------------'

echo
echo Configure complete, now type \'make\' to build.
echo

])
