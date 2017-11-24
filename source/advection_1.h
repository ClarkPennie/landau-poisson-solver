/* This is the header file associated to advection_1.cpp in which the prototypes for the functions
 * contained in that file are declared and any variables defined there to be used throughout the
 * other files are declared as external.  Any other header files which must be linked to for the
 * functions here are also included.
 *
 */

#include "LP_ompi.h"																					// allows the libraries included, macros defined and external variables declared in LP_ompi.h to be used in the advection_1 functions

//************************//
//   EXTERNAL VARIABLES   //
//************************//

extern double wt[5];																					// weights for Gaussian quadrature
extern double vt[5];																					// node values for Gaussian quadrature over the interval [-1,1]

//************************//
//   FUNCTION PROTOTYPES  //
//************************//

double Gridv(double m);

double Gridx(double m);

double rho_x(double x, double *U, int i);

double rho(double *U, int i);

double computePhi_x_0(double *U);

double computePhi(double *U, double x, int ix);

void PrintPhiVals(double *U, FILE *phifile);

double computeC_rho(double *U, int i);

double Int_Int_rho(double *U, int i);

double Int_Int_rho1st(double *U, int i);

double Int_E(double *U, int i);

double Int_E1st(double *U, int i);

double Int_fE(double *U, int i, int j);

double Int_E2nd(double *U, int i);

double I1(double *U, int k, int l);

double I2(double *U, int k, int l);

double I3(double *U, int k, int l);

double I5(double *U, int k, int l);

void computeH(double *U);

void RK3(double *U);
