/* This is the header file associated to conservationRoutines.cpp in which the prototypes for the
 * functions contained in that file are declared.  Any other header files which must be linked to
 * for the functions here are also included.
 *
 */

//************************//
//     INCLUDE GUARDS     //
//************************//

#ifndef CONSERVATIONROUTINES_H_
#define CONSERVATIONROUTINES_H_

//************************//
//        INCLUDES        //
//************************//

#include "LP_ompi.h"																					// allows the libraries included, macros defined and external variables declared in LP_ompi.h to be used in the conservationRoutines functions

//************************//
//   FUNCTION PROTOTYPES  //
//************************//

double sinc(double x);

void solveWithCCt(int nElem, double *b);

void conserveMoments(fftw_complex *qHat, fftw_complex *qHat_linear = NULL);

void createCCtAndPivot();

void conserveAllMoments(fftw_complex *qHat, fftw_complex *qHat_linear);

void conserveAllMoments_FandL(fftw_complex *qHat, fftw_complex *qHat_linear);

void conserveAllMoments_Normal(fftw_complex *qHat);

void createCCtAndPivot_AllMoments();

void conserveMass(fftw_complex *qHat, fftw_complex *qHat_linear);

void conserveMass_FandL(fftw_complex *qHat, fftw_complex *qHat_linear);

void conserveMass_Normal(fftw_complex *qHat);

void createCCtAndPivot_OnlyMass();

#endif /* CONSERVATIONROUTINES_H_ */
