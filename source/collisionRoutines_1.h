/* This is the header file associated to collisionRoutines_1.cpp in which the prototypes for the
 * functions contained in that file are declared.  Any other header files which must be linked to
 * for the functions here are also included.
 *
 */

#include "LP_ompi.h"																					// allows the libraries included, macros defined and external variables declared in LP_ompi.h to be used in the collisionRoutines_1 functions
#include "advection_1.h"																				// allows the external variables and function prototypes declared in advection_1.h to be used in the collitionRoutines_1 functions
#include "conservationRoutines.h"																		// allows the function prototypes declared in conservationRoutines.h to be used in the collisionRoutines_1 functions

//extern double IntM[10];																				// declare an array IntM to hold 10 double variables
//#pragma omp threadprivate(IntM)																	// start the OpenMP parallel construct to start the threads which will run in parallel, passing IntM to each thread as private variables which will have their contents deleted when the threads finish (doesn't seem to be doing anything since no {} afterwards???)


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
