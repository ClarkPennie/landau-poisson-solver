/* This is the header file associated to MarginalCreaetion.cpp in which the prototypes for the
 * functions contained in that file are declared.  Any other header files which must be linked to
 * for the functions here are also included.
 *
 */

//************************//
//     INCLUDE GUARDS     //
//************************//

#ifndef MARGINALCREATION_H_
#define MARGINALCREATION_H_

//************************//
//        INCLUDES        //
//************************//

#include "LP_ompi.h"																					// allows the libraries included, macros defined and external variables declared in LP_ompi.h to be used in the MarginalCreation functions
#include "advection_1.h"																				// allows the external variables and function prototypes declared in advection_1.h to be used in the MarginalCreation functions

//************************//
//   FUNCTION PROTOTYPES  //
//************************//

double f_marg(double *U, int i, int j1, double x, double v1);

void PrintMarginalLoc(FILE *margfile);

void PrintMarginal(double *U, FILE *margfile);

#endif /* MARGINALCREATION_H_ */
