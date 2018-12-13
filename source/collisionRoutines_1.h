/* This is the header file associated to collisionRoutines_1.cpp in which the prototypes for the
 * functions contained in that file are declared.  Any other header files which must be linked to
 * for the functions here are also included.
 *
 */

//************************//
//     INCLUDE GUARDS     //
//************************//

#ifndef COLLISIONROUTINES_H_
#define COLLISIONROUTINES_H_

//************************//
//        INCLUDES        //
//************************//

#include "LP_ompi.h"																					// allows the libraries included, macros defined and external variables declared in LP_ompi.h to be used in the collisionRoutines_1 functions
#include "advection_1.h"																				// allows the external variables and function prototypes declared in advection_1.h to be used in the collitionRoutines_1 functions
#include "conservationRoutines.h"																		// allows the function prototypes declared in conservationRoutines.h to be used in the collisionRoutines_1 functions

//************************//
//   FUNCTION PROTOTYPES  //
//************************//

double S1hat(double ki1,double ki2,double ki3);

double S233hat(double ki1, double ki2, double ki3);

double S213hat(double ki1, double ki2, double ki3);

double gHat3(double eta1, double eta2, double eta3, double ki1, double ki2, double ki3 );

double gHat3_linear(double eta1, double eta2, double eta3, double ki1, double ki2, double ki3 );

void generate_conv_weights(double **conv_weights);

void generate_conv_weights_linear(double **conv_weights_linear);

void fft3D(fftw_complex *in, fftw_complex *out);

void ifft3D(fftw_complex *in, fftw_complex *out);

void FS(fftw_complex *in, fftw_complex *out);

#ifdef MPI_parallelcollision
void ComputeQ(double *f, fftw_complex *qHat, double **conv_weights);
#endif

void IntModes(int k1, int k2,  int k3, int j1, int j2, int j3, double *result);

void ProjectedNodeValue(fftw_complex *qHat, double *Q_incremental);

void ComputeQ_FandL(double *f, fftw_complex *qHat, double **conv_weights, fftw_complex *qHat_linear, double **conv_weights_linear);

void RK4_FandL(double *f, int l, fftw_complex *qHat, double **conv_weights, fftw_complex *qHat_linear, double **conv_weights_linear, double *U, double *dU = NULL);

void ComputeQ(double *f, fftw_complex *qHat, double **conv_weights);

void RK4(double *f, int l, fftw_complex *qHat, double **conv_weights, double *U, double *dU = NULL);

/*
#ifdef UseMPI
	#ifdef FullandLinear
	void ComputeQ(double *f, fftw_complex *qHat, double **conv_weights, fftw_complex *qHat_linear, double **conv_weights_linear);

	void RK4(double *f, int l, fftw_complex *qHat, double **conv_weights, fftw_complex *qHat_linear, double **conv_weights_linear, double *U, double *dU);
	#else
	void ComputeQ(double *f, fftw_complex *qHat, double **conv_weights);

	void RK4(double *f, int l, fftw_complex *qHat, double **conv_weights, double *U, double *dU);
	#endif
#else
	#ifdef FullandLinear
	void ComputeQ(double *f, fftw_complex *qHat, double **conv_weights, double **conv_weights_linear);
	#else
	void ComputeQ(double *f, fftw_complex *qHat, double **conv_weights);

	void RK4(double *f, int l, fftw_complex *qHat, double **conv_weights, double *U);
	#endif
#endif  
*/
#endif /* COLLISIONROUTINES_H_ */
