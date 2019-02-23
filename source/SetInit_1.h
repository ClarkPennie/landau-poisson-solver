/* This is the header file associated to SetInit_1.cpp in which the prototypes for the
 * functions contained in that file are declared.  Any other header files which must be linked to
 * for the functions here are also included.
 *
 */

//************************//
//     INCLUDE GUARDS     //
//************************//

#ifndef SETINIT_1_H_
#define SETINIT_1_H_

//************************//
//        INCLUDES        //
//************************//

#include "LP_ompi.h"																					// allows the libraries included, macros defined and external variables declared in LP_ompi.h to be used in the SetInit_1 functions
#include "advection_1.h"																				// allows the external variables and function prototypes declared in advection_1.h to be used in the SetInit_1 functions
#include "FieldCalculations.h"																			// allows the function prototypes declared in FieldCalculations.h to be used in the advection_1 functions

//************************//
//   FUNCTION PROTOTYPES  //
//************************//

void trapezoidalRule(int nPoints, double *weight);

double f_TS(double v1, double v2, double v3);

double f_2Gauss(double v1, double v2, double v3);

double Mw(double v1, double v2, double v3, double T);

double Mw_x(double x);

double f_2H(double x);

void SetInit_LD(double *U);

void SetInit_4H(double *U, double T0, double C);

void SetInit_2H(double *U);

void SetInit_ND(double *U);

void setInit_spectral(double *U, double **f);

#endif /* SETINIT_1_H_ */
