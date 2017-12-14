/* This is the header file associated to NegativityChecks.cpp in which the prototypes for the
 * functions contained in that file are declared.  Any other header files which must be linked to
 * for the functions here are also included.
 *
 */

//************************//
//     INCLUDE GUARDS     //
//************************//

#ifndef NEGATIVITYCHECKS_H_
#define NEGATIVITYCHECKS_H_

//************************//
//        INCLUDES        //
//************************//

#include "LP_ompi.h"																					// allows the libraries included, macros defined and external variables declared in LP_ompi.h to be used in the NegativityChecks functions
#include "advection_1.h"																				// allows the external variables and function prototypes declared in advection_1.h to be used in the NegativityChecks functions

//************************//
//   FUNCTION PROTOTYPES  //
//************************//

double computeCellAvg(double *U, int i, int j1, int j2, int j3);

void FindNegVals(double *U, int *NegVals, double *AvgVals);

void FindNegVals_old(double *U, int *NegVals);

void CheckNegVals(double *U, int *NegVals, double *AvgVals);

#endif /* NEGATIVITYCHECKS_H_ */
