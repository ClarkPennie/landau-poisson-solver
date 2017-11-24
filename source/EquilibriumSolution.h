/* This is the header file associated to EquilibriumSolutions.cpp in which the prototypes for the
 * functions contained in that file are declared.  Any other header files which must be linked to
 * for the functions here are also included.
 *
 */

#include "LP_ompi.h"																					// allows the libraries included, macros defined and external variables declared in LP_ompi.h to be used in the EquilibriumSolution functions
#include "advection_1.h"																				// allows the external variables and function prototypes declared in advection_1.h to be used in the EquilibriumSolution functions

//************************//
//   FUNCTION PROTOTYPES  //
//************************//

void ExportRhoQuadVals(double *U);

void ComputeEquiVals(double *f);

void PrintEquiVals(double *U, FILE *margfile);
