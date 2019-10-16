/* This is the header file associated to FieldCalculations.cpp in which the prototypes for the functions
 * contained in that file are declared.  Any other header files which must be linked to for the
 * functions here are also included.
 *
 */

//************************//
//     INCLUDE GUARDS     //
//************************//

#ifndef FIELDCALCULATIONS_H_
#define FIELDCALCULATIONS_H_

//************************//
//        INCLUDES        //
//************************//

#include "LP_ompi.h"																					// allows the libraries included, macros defined and external variables declared in LP_ompi.h to be used in the advection_1 functions
#include "advection_1.h"																				// allows the external variables and function prototypes declared in advection_1.h to be used in the FieldCalculations functions


//************************//
//   FUNCTION PROTOTYPES  //
//************************//

double rho_x(double x, double *U, int i);

double rho(double *U, int i);

void PrintFieldLoc(FILE *phifile, FILE *Efile);

void PrintFieldData(double* U, FILE *phifile, FILE *Efile);

double computeC_rho(double *U, int i);

double Int_Int_rho(double *U, int i);

double Int_Int_rho1st(double *U, int i);

double Int_fE(double *U, int i, int j);

double computePhi_x_0(double *U);

double computePhi(double *U, double x, int ix);

double computeE(double *U, double x, int ix);

double Int_E(double *U, int i);

double Int_E1st(double *U, int i);

double Int_E2nd(double *U, int i);

void PrintFieldData(double* U_vals, FILE *phifile, FILE *Efile);

double Int_E(double *U, int i);

double Int_E1st(double *U, int i);

double Int_E2nd(double *U, int i);

double computePhi_x_0_Normal(double *U);

double computePhi_Normal(double *U, double x, int ix);

void PrintFieldData_Normal(double* U_vals, FILE *phifile, FILE *Efile);

double Int_E_Normal(double *U, int i);

double Int_E1st_Normal(double *U, int i);

double Int_E2nd_Normal(double *U, int i);

double DopingProfile(int i);

void PrintDoping();

double EpsilonValue(int i);

void PrintEpsilon();

double computePhi_x_0_Doping(double *U);

double computePhi_x_0_Doping(double *U);

double computePhi_Doping(double *U, double x, int ix);

double computeE_Doping(double *U, double x, int ix);

void PrintFieldData_Doping(double* U_vals, FILE *phifile, FILE *Efile);

double Int_E_Doping(double *U, int i);

double Int_E1st_Doping(double *U, int i);

double Int_E2nd_Doping(double *U, int i);

#endif /* FIELDCALCULATIONS_H_ */

