/* This is the header file associated to advection_1.cpp in which the prototypes for the functions
 * contained in that file are declared and any variables defined there to be used throughout the
 * other files are declared as external.  Any other header files which must be linked to for the
 * functions here are also included.
 *
 */

//************************//
//     INCLUDE GUARDS     //
//************************//

#ifndef ADVECTION_1_H_
#define ADVECTION_1_H_

//************************//
//        INCLUDES        //
//************************//

#include "LP_ompi.h"																					// allows the libraries included, macros defined and external variables declared in LP_ompi.h to be used in the advection_1 functions
#include "FieldCalculations.h"																			// allows the function prototypes declared in FieldCalculations.h to be used in the advection_1 functions
#include "SetInit_1.h"																					// allows the function prototypes declared in SetInit_1.h to be used in the advection_1 functions

//************************//
//   EXTERNAL VARIABLES   //
//************************//

extern double wt[5];																					// weights for Gaussian quadrature
extern double vt[5];																					// node values for Gaussian quadrature over the interval [-1,1]

//************************//
//   FUNCTION PROTOTYPES  //
//************************//

double Gridv(double m);

double Gridv_L(double m);

double Gridv_H(double m);

double Gridx(double m);

void DirichletBC(vector<double>& Ub_vals, int i, int j1, int j2, int j3);

double I1(double *U, int k, int l);

double I2(double *U, int k, int l);

double I3_B(double *U, int k, int l);

double I3(double *U, int k, int l);

double I3_Doping(double *U, int k, int l);

double I3_Normal(double *U, int k, int l);

double I5(double *U, int k, int l);

void computeH(double *U);

void RK3(double *U);

#endif /* ADVECTION_1_H_ */
