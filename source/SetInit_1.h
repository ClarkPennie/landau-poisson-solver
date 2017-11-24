/* This is the header file associated to SetInit_1.cpp in which the prototypes for the
 * functions contained in that file are declared.  Any other header files which must be linked to
 * for the functions here are also included.
 *
 */

#include "LP_ompi.h"																					// allows the libraries included, macros defined and external variables declared in LP_ompi.h to be used in the SetInit_1 functions
#include "advection_1.h"																				// allows the external variables and function prototypes declared in advection_1.h to be used in the SetInit_1 functions

//************************//
//   FUNCTION PROTOTYPES  //
//************************//

void trapezoidalRule(int nPoints, double *weight);

double f_TS(double v1, double v2, double v3);

double f_2Gauss(double v1, double v2, double v3);

double Mw(double v1, double v2, double v3);

double Mw_x(double x);

double f_2H(double x);

void SetInit_LD(double *U);

void SetInit_4H(double *U);

void SetInit_2H(double *U);

#ifdef MPI
void setInit_spectral(double *U, double **f);
#else
void setInit_spectral(double *U, double **f);
#endif
