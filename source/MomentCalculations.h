/* This is the header file associated to MomentCalculations.cpp in which the prototypes for the
 * functions contained in that file are declared.  Any other header files which must be linked to
 * for the functions here are also included.
 *
 */

//************************//
//     INCLUDE GUARDS     //
//************************//

#ifndef MOMENTCALCULATIONS_H_
#define MOMENTCALCULATIONS_H_

//************************//
//        INCLUDES        //
//************************//

#include "LP_ompi.h"																					// allows the libraries included, macros defined and external variables declared in LP_ompi.h to be used in the MomentCalculations functions
//#include "advection_1.h"																				// allows the external variables and function prototypes declared in advection_1.h to be used in the MomentCalculations functions
#include "FieldCalculations.h"																			// allows the function prototypes declared in FieldCalculations.h to be used in the advection_1 functions

//************************//
//   FUNCTION PROTOTYPES  //
//************************//

double computeMass(double *U);

double computeMass_Inhomo(double *U);

double computeMass_Inhomo_Uniform(double *U);

double computeMass_Inhomo_Refined(double *U);

double computeMass_Homo(double *U);

double computeMass_in_x(double *U, int i, double x);

void computeMomentum(double *U, double *a);

void computeMomentum_Inhomo(double *U, double *a);

void computeMomentum_Inhomo_Uniform(double *U, double *a);

void computeMomentum_Inhomo_Refined(double *U, double *a);

void computeMomentum_Homo(double *U, double *a);

double computeBulkVelocity_v1_in_x(double *U, int i, double x);

double computeBulkVelocity_v2_in_x(double *U, int i, double x);

double computeBulkVelocity_v3_in_x(double *U, int i, double x);

double computeBulkMomentum_nu1_in_x(double *U, int i, double x);

double computeBulkMomentum_nu2_in_x(double *U, int i, double x);

double computeBulkMomentum_nu3_in_x(double *U, int i, double x);

double computeKiE(double *U);

double computeKiE_Inhomo_Uniform(double *U);

double computeKiE_Inhomo_Refined(double *U);

double computeKiE_Inhomo(double *U);

double computeKiE_Homo(double *U);

double computeT_in_x(double *U, int i, double x);

double computeKiE_in_x(double *U, int i, double x);

double computeKiEratio(double *U, int *NegVals);

double computeKiEratio_Inhomo(double *U, int *NegVals);

double computeKiEratio_Homo(double *U, int *NegVals);

double computeEleE(double *U);

#endif /* MOMENTCALCULATIONS_H_ */
