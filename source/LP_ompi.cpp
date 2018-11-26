/* This is the main source file for running the simulations of the space inhomogeneous Landau-Poisson
 * equation, using the conservative spectral method with RK4 for the collision problem and RKDG for
 * the advection problem (both of which result from time-splitting).
 *
 */

//***********************************************//
//                                               //
//          MACROS LOCATED IN LP_ompi.h          //
//                                               //
//      Use these to choose:                     //
//      - whether or not to use MPI,             //
//      - the initial conditions,                //
//      - if this is the first run or not.       //
//                                               //
//***********************************************//

#include "LP_ompi.h"																				// LP_ompi.h is where the libraries required by the program included, all macros (to decide the behaviour of a given run) are defined and all variables to be used throughout the various files are defined as external

double PI=M_PI;																						// declare PI and set it to M_PI (the value stored in the library math.h)
int M=5;																							// declare M (the number of collision invarients) and set it equal to 5
int Nx, Nv, nT, N;			 											 							// declare Nx (no. of x discretised points), Nv (no. of v discretised point), nT (no. of time discretised points) & N (no. of nodes in the spectral method) and setting all their values
int size_v, size, size_ft; 																			// declare size_v (no. of total v discretised points in 3D), size (the total no. of discretised points) & size_ft (total no. of spectral discretised points in 3D)
double dv, dx; 																						// declare dv (the velocity stepsize) & dx (the space stepsize)
double A_amp, k_wave;																				// declare A_amp (the amplitude of the perturbing wave) & k_wave (the wave number of the perturbing wave)
double Lx, Lv;																						// declare Lx (for 0 < x < Lx) and set it to & Lv (for -Lv < v < Lv in the advection problem)
double L_v, R_v, L_eta;																				// declare L_v (for -Lv < v < Lv in the collision problem), R_v (for v in B_(R_v) in the collision problem) & L_eta (for Fourier space, -L_eta < eta < L_eta)
double h_eta, h_v;																					// declare h_eta (the Fourier stepsize) & h_v (also the velocity stepsize but for the collision problem)
double nu, dt; 																						// declare nu (1/knudson#), dt (the timestep) & nthread (the number of OpenMP threads)

double *v, *eta;																					// declare v (the velocity variable) & eta (the Fourier space variable)
double *wtN;																						// declare wtN (the trapezoidal rule weights to be used)
double scale, scale3, scaleL, scalev;																// declare scale (the 1/sqrt(2pi) factor appearing in Gaussians), scale (the 1/(sqrt(2pi))^3 factor appearing in the Maxwellian), scaleL (the volume of the velocity domain) & scalev (the volume of a discretised velocity element)
double **C1, **C2;																					// declare pointers to matrices C1 (the real part of the conservation matrix C) & C2 (the imaginary part of the conservation matrix C), CCt (of dimension 5x5) & CCt_linear (of dimension 2x2)
double CCt[5*5], CCt_linear[2*2];																	// declare matrices CCt (C*C^T, for the conservation matrix C) & CCt_linear (C*C^T, for the conservation matrix C, in the two species collision operator)
double lamb[5], lamb_linear[2];																		// declare the arrays lamb (to hold 5 values) & lamb_linear (to hold 2 values)

double *U1, *Utmp, *output_buffer_vp;//, **H;														// declare pointers to U1, Utmp (both used to help store the values in U, declared later) & output_buffer_vp (a buffer used when sending the data between MPI processes during the VP method)
double *Q, *f1, *Q1, *Utmp_coll;//*f2, *f3;															// declare pointers to Q (the discretised collision operator), f1 (used to help store the solution during the collisional problem), Q1 (used in calculation of the collision operator) & Utmp_coll (used to store calculations from the RK4 method used in the collisional problem)
fftw_complex *Q1_fft, *Q2_fft, *Q3_fft;																// declare pointers to the complex numbers Q1_fft, Q2_fft & Q3_fft (involved in storing the FFT of Q)

#ifdef FullandLinear																				// only do this if FullandLinear was defined
fftw_complex *Q1_fft_linear, *Q2_fft_linear, *Q3_fft_linear;										// declare pointers to the complex numbers Q1_fft_linear, Q2_fft_linear & Q3_fft_linear (involved in storing the FFT of the two species collison operator Q)
#endif

fftw_complex *fftIn, *fftOut;																		// declare pointers to the FFT variables fftIn (a vector to be to have the FFT applied to it) & fftOut (the output of an FFT)

double ce, *cp, *intE, *intE1, *intE2;																// declare ce and pointers to cp, intE, intE1 & intE2 (precomputed quantities for advections)

// SET UP FFT PLANS (WHICH ARE USED MULTIPLE TIMES):
fftw_plan p_forward; 																				// declare the fftw_plan p_forward (an object which contains all the data which allows fftw3 to compute the FFT)
fftw_plan p_backward; 																				// declare the fftw_plan p_backward (an object which contains all the data which allows fftw3 to compute the inverse FFT)
fftw_complex *temp;																					// declare a pointer to complex number temp (a temporary array used when calculating the FFT)

int myrank_mpi, nprocs_mpi, nprocs_Nx;																// declare myrank_mpi (the rank of the current MPI process running), nprocs_mpi (the total number of MPI processes) & nprocs_Nx (the amount of MPI processes used for the collisionless VP problem)
int chunksize_dg, chunksize_ft, chunk_Nx;															// declare chunksize_dg (the amount of data each processes works on during the DG method in the VP problem), chunksize_ft (the amount of data each process works on during the collisional problem) & chunk_Nx (the number of space-steps sent to each process in the collisional problem)

int *fNegVals;																						// declare fNegVals (to store where DG solution goes negative - a 1 if negative and a 0 if positive)
double *fAvgVals;																					// declare fAvgVals (to store the average values of f on each cell)
double *fEquiVals;																					// declare f_equivals (to store the equilibrium solution)

int main()
{
	int i, j, k, j1, j2, j3, l; 																	// declare i, j, k (counters), j1, j2, j3 (velocity space counters) & l (the index of the current DG basis function being integrated against)
	int  tp, t=0; 																					// declare tp (the amount size of the data which stores the DG coefficients of the solution read from a previous run) & t (the current time-step) and set it to 0
	int k_v, k_eta, k_local, nprocs_vlasov;															// declare k_v (the index of a DG coefficient), k_eta (the index of a DG coefficient in Fourier space), k_local (the index of a DG coefficient, local to the space chunk on the current process) & nprocs_vlasov (the number of processes used for solving the Vlasov equation)
	double tmp, mass, a[3], KiE, EleE, KiEratio, ent1, l_ent1, ll_ent1;								// declare tmp (the square root of electric energy), mass (the mass/density rho), a (the momentum vector J), KiE (the kinetic energy), EleE (the electric energy), KiEratio (the ratio of kinetic energy between where f is positive and negative),  ent1 (the entropy with negatives discarded), l_ent1 (log of the ent1) & ll_ent1 (log of log of ent1)
	double *U, **f, *output_buffer;//, **conv_weights_local;										// declare pointers to U (the vector containing the coefficients of the DG basis functions for the solution f(x,v,t) at the given time t), f (the solution which has been transformed from the DG discretisation to the appropriate spectral discretisation) & output_buffer (from where to send MPI messages)
	double **conv_weights, **conv_weights_linear;													// declare a pointer to conv_weights (a matrix of the weights for the convolution in Fourier space of single species collisions) conv_weights_linear (a matrix of convolution weights in Fourier space of two species collisions)
  
	fftw_complex *qHat, *qHat_linear;																// declare pointers to the complex numbers qHat (the DFT of Q) & qHat_linear (the DFT of the two species colission operator Q);
  
	GRVY_Input_Class iparse;																		// Define a GRVY input parsing object

	if(! iparse.Open("./LPsolver-input.txt"))														// Initialise and read in the GRVY file with the input paramters
	{
		exit(1);																					// Exit if it does not exist
	}
	if (! iparse.Read_Var("nT",&nT) )																// Check if the variable nT has been set and, if not, exit
	{
		exit(1);
	}
	printf("--> %-11s = %d\n","nT",nT);																// Print the set value
	if (! iparse.Read_Var("Nx",&Nx) )																// Check if the variable Nx has been set and, if not, exit
	{
		exit(1);
	}
	printf("--> %-11s = %d\n","Nx",Nx);																// Print the set value
	if (! iparse.Read_Var("Nv",&Nv) )																// Check if the variable Nv has been set and, if not, exit
	{
		exit(1);
	}
	printf("--> %-11s = %d\n","Nv",Nv);																// Print the set value
	if (! iparse.Read_Var("N",&N) )																	// Check if the variable N has been set and, if not, exit
	{
		exit(1);
	}
	printf("--> %-11s = %d\n","N",N);																// Print the set value
	if (! iparse.Read_Var("nu",&nu) )																// Check if the variable dt has been set and, if not, exit
	{
		exit(1);
	}
	printf("--> %-11s = %g\n","nu",nu);																// Print the set value
	if (! iparse.Read_Var("dt",&dt) )																// Check if the variable dt has been set and, if not, exit
	{
		exit(1);
	}
	printf("--> %-11s = %g\n","dt",dt);																// Print the set value

	#ifdef Damping																					// only do this if Damping was defined
	if (! iparse.Read_Var("Damping/A_amp",&A_amp) )													// Check if the variable A_amp has been set and, if not, exit
	{
		exit(1);
	}
	if (! iparse.Read_Var("Damping/k_wave",&k_wave) )												// Check if the variable k_wave has been set and, if not, exit
	{
		exit(1);
	}
	if (! iparse.Read_Var("Damping/Lv",&Lv) )														// Check if the variable Lv has been set and, if not, exit
	{
		exit(1);
	}
	iparse.Read_Var("Damping/Lx",&Lx,2*PI/k_wave);													// Check if the variable Lx has been set and, if not, set to the default value of 2pi/k_wave
	#endif

	#ifdef TwoStream																				// only do this if Damping was defined
	if (! iparse.Read_Var("TwoStream/A_amp",&A_amp) )												// Check if the variable A_amp has been set and, if not, exit
	{
		exit(1);
	}
	if (! iparse.Read_Var("TwoStream/Lv",&Lv) )														// Check if the variable Lv has been set and, if not, exit
	{
		exit(1);
	}
	iparse.Read_Var("TwoStream/Lx",&Lx,2*PI/k_wave);													// Check if the variable Lx has been set and, if not, set to the default value of 2pi/k_wave
	k_wave=2*PI/4.																					// set k_wave to pi/2
	#endif

	#ifdef FourHump																					// only do this if Damping was defined
	iparse.Read_Var("FourHump/A_amp",&A_amp,0.);														// Check if the variable A_amp has been set and, if not, set the value to 0 (it's not needed)
	iparse.Read_Var("FourHump/k_wave",&k_wave,0.5);													// Check if the variable k_wave has been set and, if not, set to the default value of 0.5 (to calculate Lx)
	if (! iparse.Read_Var("FourHump/Lv",&Lv) )														// Check if the variable Lv has been set and, if not, exit
	{
		exit(1);
	}
	iparse.Read_Var("FourHump/Lx",&Lx,2*PI/k_wave);													// Check if the variable Lx has been set and, if not, set to the default value of 2pi/k_wave
	#endif

	#ifdef TwoHump																					// only do this if Damping was defined
	iparse.Read_Var("TwoHump/A_amp",&A_amp,0.);														// Check if the variable A_amp has been set and, if not, set the value to 0 (it's not needed)
	iparse.Read_Var("TwoHump/k_wave",&k_wave,0.5);													// Check if the variable k_wave has been set and, if not, set to the default value of 0.5 (to calculate Lx)
	if (! iparse.Read_Var("TwoHump/Lv",&Lv) )														// Check if the variable Lv has been set and, if not, exit
	{
		exit(1);
	}
	iparse.Read_Var("TwoHump/Lx",&Lx,2*PI/k_wave);													// Check if the variable Lx has been set and, if not, set to the default value of 2pi/k_wave
	#endif

	printf("--> %-11s = %g\n","A_amp",A_amp);														// Print the value set for A_amp
	printf("--> %-11s = %g\n","k_wave",k_wave);														// Print the value set for k_wave
	printf("--> %-11s = %g\n","Lv",Lv);																// Print the value set for Lv
	printf("--> %-11s = %g\n","Lx",Lx);																// Print the value set for Lx

	size_v=Nv*Nv*Nv;																				// set size_v to Nv^3
	size=Nx*size_v;																					// set size to size_v*Nx
	size_ft=N*N*N; 																					// set size_ft to N^3
	dv=2.*Lv/Nv;																					// set dv to 2Lv/Nv
	dx=Lx/Nx; 																						// set dx to Lx/Nx
	scalev=dv*dv*dv;																				// set scalev to dv^3
	L_v=Lv;																							// set L_v to Lv
	R_v=Lv;																							// set R_v to Lv
	scaleL=8*Lv*Lv*Lv;																				// set scaleL to 8Lv^3

	//************************
	//MPI-related variables!
	//************************
	int required=MPI_THREAD_MULTIPLE;																// declare required and set it to MPI_THREAD_MULTIPLE (so that in the hybrid OpenMP/MPI routines, multiple threads may call MPI, with no restrictions)           . MPI_THREAD_SERIALIZED; // Required level of MPI threading support
	int provided;                       															// declare provided (the actual provided level of MPI thread support)
	MPI_Status status;																				// declare status (to store any status required for MPI operations)

	MPI_Init_thread(NULL, NULL, required, &provided);												// initialise the hybrid MPI & OpenMP environment, requesting the level of thread support to be required and store the actual thread support provided in provided
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);														// store the rank of the current process in the MPI_COMM_WORLD communicator in myrank_mpi
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);														// store the total number of processes running in the MPI_COMM_WORLD communicator in nprocs_mpi

	int nthread;																					// declare nthread (the number of OpenMP threads)
	nthread = omp_get_max_threads();																// set nthread to the value of the environment variable OMP_NUM_THREADS by calling the OpenMP function omp_get_max_threads
  
	// CHECK THE LEVEL OF THREAD SUPPORT:
	if (provided < required)																		// only do this if the required thread support was not possible
	{
		if (myrank_mpi == 0)
		{
			printf( "provided=%d < required=%d. Warning:  This MPI implementation "
						"provides insufficient threading support.\n", provided, required);			// have the process with rank 0 display in the output file that the required thread support was not possible
			MPI_Finalize();																			// ensure that MPI exits cleanly
			exit(0);
		}
	}
	else
	{
		omp_set_num_threads(nthread);																// if the thread support required is possible, set the number of OpenMP threads to nthread
	}
   
	double MPIt1, MPIt2, MPIelapsed;																// declare MPIt1 (the start time of an MPI operation), MPIt2 (the end time of an MPI operation) and MPIelapsed (the total time for the MPI operation)

	if(size_v%nprocs_mpi != 0)																		// check that size_v/nprocs_mpi has no remainder
	{
		if (myrank_mpi==0)
		{
			printf("Error: There must be integer number of chunks\n");								// if nprocs_mpi isn't a factor of size_v then display this in the output file
		}
		MPI_Finalize();																				// ensure that MPI exits cleanly
		exit(0);
	}

	chunksize_dg = size/nprocs_mpi;																	// set chunksize_dg to size/nprocs_mpi (which will be an integer since nprocs_mpi is a facot of size_v and size = Nx*size_v)
	chunksize_ft = size_ft/nprocs_mpi; 																// set chunksize_ft to size_ft/nprocs_mpi

	if(Nx%nprocs_mpi == 0)
	{
		chunk_Nx = Nx/nprocs_mpi;																	// if nprocs_mpi divides into Nx, set chunk_Nx to Nx/nprocs_mpi
	}
	else
	{
		chunk_Nx = Nx/nprocs_mpi + 1;																// if nprocs_mpi does not divide into Nx, set chunk_Nx to Nx/nprocs_mpi + 1
	}

	nprocs_Nx = (int)((double)Nx/(double)chunk_Nx + 0.5);											// set nprocs_Nx to Nx/chunk_Nx + 0.5 and store the result as an integer
	
	U = (double*)malloc(size*6*sizeof(double));														// allocate enough space at the pointer U for 6*size many double numbers
	U1 = (double*)malloc(size*6*sizeof(double));													// allocate enough space at the pointer U1 for 6*size many floating point numbers
 
	Utmp = (double*)malloc(chunksize_dg*6*sizeof(double));											// allocate enough space at the pointer Utmp for 6*chunksize_dg many floating point numbers
	// H[i] = (double*)malloc(6*sizeof(double));}
	output_buffer_vp = (double *) malloc(chunksize_dg*6*sizeof(double));							// allocate enough space at the pointer output_buffer_vp for 6*chunksize_dg many floating point numbers
  
	cp = (double*)malloc(Nx*sizeof(double));														// allocate enough space at the pointer cp for Nx many double numbers
	intE = (double*)malloc(Nx*sizeof(double));														// allocate enough space at the pointer intE for Nx many double numbers
	intE1 = (double*)malloc(Nx*sizeof(double));														// allocate enough space at the pointer intE1 for Nx many double numbers
	intE2 = (double*)malloc(Nx*sizeof(double));														// allocate enough space at the pointer intE2 for Nx many double numbers

	fNegVals = (int*)malloc(size*sizeof(int));														// allocate enough space at the pointer fNegVals for size many integers
	fAvgVals = (double*)malloc(size*sizeof(double));												// allocate enough space at the pointer fAvgVals for size many doubles
	fEquiVals = (double*)malloc(5*5*5*5*size*sizeof(double));										// allocate enough space at the pointer fEquiVals for 5*Nx*(5*Nv)^3 many doubles

	if(nu > 0.)
	{
		C1 = (double**)malloc(M*sizeof(double *)); 													// allocate enough space at the pointer C1 for M many pointers to double numbers
		C2 = (double**)malloc(M*sizeof(double *));													// allocate enough space at the pointer C2 for M many pointers to double numbers
		for(i=0;i<M;i++)
		{
			C1[i] = (double *)malloc(size_ft*sizeof(double));										// allocate enough space at the ith entry of C1 for size_ft many double numbers
			C2[i] = (double *)malloc(size_ft*sizeof(double));										// allocate enough space at the ith entry of C2 for size_ft many double numbers
		}
		f = (double **)malloc(chunk_Nx*sizeof(double *));											// allocate enough space at the pointer f for chunk_Nx many pointers to double numbers
		for (i=0;i<chunk_Nx;i++)
		{
			f[i] = (double *)malloc(size_ft*sizeof(double));										// allocate enough space at the ith entry of f for size_ft many double numbers
		}
		conv_weights = (double **)malloc(size_ft*sizeof(double *));									// allocate enough space at the pointer conv_weights for size_ft many pointers to float numbers
		for (i=0;i<size_ft;i++)
		{
			conv_weights[i] = (double *)malloc(size_ft*sizeof(double));								// allocate enough space at the ith entry of conv_weights for size_ft many float numbers
		}

		#ifdef FullandLinear																		// only do this if FullandLinear was defined
		conv_weights_linear = (double **)malloc(size_ft*sizeof(double *));							// allocate enough space at the pointer conv_weight_linear for size_ft many pointers to float numbers
		for (i=0;i<size_ft;i++)
		{
			conv_weights_linear[i] = (double *)malloc(size_ft*sizeof(double));						// allocate enough space at the ith entry of conv_weights_linear for size_ft many float numbers
		}
  
		Q1_fft_linear = (fftw_complex *)fftw_malloc(size_ft*sizeof(fftw_complex));					// allocate enough space at the pointer Q1_fft_linear for size_ft many complex numbers
		Q2_fft_linear = (fftw_complex *)fftw_malloc(size_ft*sizeof(fftw_complex));					// allocate enough space at the pointer Q2_fft_linear for size_ft many complex numbers
		Q3_fft_linear = (fftw_complex *)fftw_malloc(size_ft*sizeof(fftw_complex));					// allocate enough space at the pointer Q3_fft_linear for size_ft many complex numbers
  
		qHat_linear = (fftw_complex *)fftw_malloc(size_ft*sizeof(fftw_complex));					// allocate enough space at the pointer QHat_linear for size_ft many complex numbers
		#endif

		Q = (double*)malloc(size_ft*sizeof(double));												// allocate enough space at the pointer Q for size_ft many double numbers
		f1 = (double*)malloc(size_ft*sizeof(double)); 												// allocate enough space at the pointer f1 for size_ft many double numbers
		Q1 = (double*)malloc(size_ft*sizeof(double));												// allocate enough space at the pointer Q1 for size_ft many double numbers
		Utmp_coll = (double*)malloc(chunk_Nx*size_v*5*sizeof(double));								// allocate enough space at the pointer Utmp_coll for 5*chunk_Nx*size_v many double numbers
		output_buffer = (double*)malloc(chunk_Nx*size_v*5*sizeof(double));							// allocate enough space at the pointer output_buffer for 5*chunk_Nx*size_v many double numbers
  
		//f2 = (double *)malloc(size_ft*sizeof(double));
		//Q2 = (double *)malloc(N*N*N*sizeof(double));
		//f3 = (double *)malloc(size_ft*sizeof(double));
		//Q3 = (double *)malloc(N*N*N*sizeof(double));

		temp = (fftw_complex *)fftw_malloc(size_ft*sizeof(fftw_complex));							// allocate enough space at the pointer temp for size_ft many complex numbers
		//qHat_local = (fftw_complex *)fftw_malloc(chunksize_ft*sizeof(fftw_complex));
		qHat = (fftw_complex *)fftw_malloc(size_ft*sizeof(fftw_complex));							// allocate enough space at the pointer qHat for size_ft many complex numbers
  
		// INITIALISE FFTW FOR USE WITH THREADING (MUST BE DONE BEFORE ANY PLAN IS CREATED):
		fftw_init_threads();																		// initialise the environment for using the fftw3 routines with multiple threads
		fftw_plan_with_nthreads(nthread);															// set the number of threads used by fftw3 routines to nthread

		// SET UP PLANS FOR FFTs (EXECUTED BY USING nThreads):
		p_forward = fftw_plan_dft_3d (N, N, N, temp, temp, FFTW_FORWARD, FFTW_MEASURE);				// set p_forward to a 3D fftw plan of dimension NxNxN, which will take the FFT of the vector in temp, store the result back in temp, set the sign to FFTW_FORWARD (so that this is an FFT) and set the flag to FFT_MEASURE so that at this stage fftw3 finds the most efficient way to compute the FFT of this size
		p_backward = fftw_plan_dft_3d (N, N, N, temp, temp, FFTW_BACKWARD, FFTW_MEASURE);			// set p_backward to a 3D fftw plan of dimension NxNxN, which will take the FFT of the vector in temp, store the result back in temp, set the sign to FFTW_BACKWARD (so that this is an inverse FFT) and set the flag to FFT_MEASURE so that at this stage fftw3 finds the most efficient way to compute the FFT of this size

		wtN = (double *)malloc(N*sizeof(double));													// allocate enough space at the pointer wtN to store N many double numbers
		v = (double *)malloc(N*sizeof(double));														// allocate enough space at the pointer v to store N many double numbers
		eta = (double *)malloc(N*sizeof(double));													// allocate enough space at the pointer eta to store N many double numbers
  
		Q1_fft = (fftw_complex *)fftw_malloc(size_ft*sizeof(fftw_complex));							// allocate enough space at the pointer Q1_fft for size_ft many complex numbers
		Q2_fft = (fftw_complex *)fftw_malloc(size_ft*sizeof(fftw_complex));							// allocate enough space at the pointer Q2_fft for size_ft many complex numbers
		Q3_fft = (fftw_complex *)fftw_malloc(size_ft*sizeof(fftw_complex));							// allocate enough space at the pointer Q3_fft for size_ft many complex numbers
  
		fftOut = (fftw_complex *)fftw_malloc(size_ft*sizeof(fftw_complex));							// allocate enough space at the pointer fftOut for size_ft many complex numbers
		fftIn = (fftw_complex *)fftw_malloc(size_ft*sizeof(fftw_complex));							// allocate enough space at the pointer fftIn for size_ft many complex numbers
 
  
		char buffer_weights[100], loading_buffer[100];												// declare the arrays buffer_weights (to store a string which displays the values of N & L_v) & loading_buffer (...?)
		sprintf(buffer_weights,"Weights/N%d_L%g_Landau.wts",N,L_v);									// store the values of N & L_v in buffer_weights
		#ifdef FullandLinear																		// only do this if Fullandlinear was defined
		char buffer_weights1[100];																	// declare the array buffer_weights1 (to store a string which displays the values of N & L_v for the linear case)
		sprintf(buffer_weights1,"Weights/N%d_L%g_Landau_linear.wts",N,L_v);							// store the values of N & L_v in buffer_weights1 and note that it's for the linear Landau damping
		#endif

		// COMMONLY USED CONSTANTS:
		scale = 1.0/sqrt(2.0*M_PI);																	// set scale to 1/sqrt(2*pi)
		scale3 = pow(scale, 3.0);																	// set scale3 to scale^3

		//INITIALISE VELOCITY AND FOURIER DOMAINS:
		L_v = Lv;//sqrt(0.5*(double)N*PI); 															// set L_v to Lv
		L_eta = 0.5*(double)(N-1)*PI/L_v;					// BUG: N-1?							// set L_eta to (N-1)*Pi/(2*L_v)
		h_v = 2.0*L_v/(double)(N-1);						// BUG: N-1?							// set h_v to 2*L_v/(N-1)
		h_eta = 2.0*L_eta/(double)(N);						// BUG: N?								// set h_eta to 2*L_eta/N
		for(i=0;i<N;i++)																			// store the discretised velocity and Fourier space points
		{
			eta[i] = -L_eta + (double)i*h_eta;														// set the ith value of eta to -L_eta + i*h_eta
			v[i] = -L_v + (double)i*h_v;															// set the ith value of v to -L_v + i*h_v
		}

		trapezoidalRule(N, wtN);																	// set wtN to the weights required for a trapezoidal rule with N points
  
		createCCtAndPivot();																		// calculate the values of the conservation matrices

		FILE *fidWeights;																			// declare a pointer to the file fidWeights (which will store the precomputed weights in Fourier transform of Q)
		//check to see if these convolution weights have already been pre-computed
		/*if(myrank_mpi==0){
		if(fidWeights = fopen(buffer_weights,"r")) { //checks if we've already stored a file
		printf("Stored weights found. Reading in weights... \n");
			for(i=0;i<size_ft;i++)
			{
				tp=fread(conv_weights[i], sizeof(double), size_ft, fidWeights);
				if(tp != size_ft){printf("reaing error! tp=%d\n", tp);exit(0);}
			}
		  fclose(fidWeights);
			}
			else {
		  printf("Stored weights NOT found. Please go back and launch MPI_WeightGenerator. \n");
		  MPI_Finalize();
		  exit (1);
			}
		  }

		  for(k=0;k<size_ft;k++)MPI_Bcast(conv_weights[k], size_ft, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  MPI_Barrier(MPI_COMM_WORLD);
  
		  #ifdef FullandLinear
		  FILE *fidWeights1;
		  if(myrank_mpi==0){
			if(fidWeights1 = fopen(buffer_weights1,"r")) { //checks if we've already stored a file
		  printf("Stored weights found. Reading in weights... \n");
		  for(i=0;i<size_ft;i++) {
			tp=fread(conv_weights_linear[i], sizeof(double), size_ft, fidWeights1);
			if(tp != size_ft){printf("reaing error! tp=%d\n", tp);exit(0);}
		  }
		  fclose(fidWeights1);
			}
			else {
			  printf("Stored weights NOT found. Please go back and launch MPI_WeightGenerator. \n");
			  MPI_Finalize();
			  exit (1);
			}
		  }
  
		  for(k=0;k<size_ft;k++)MPI_Bcast(conv_weights_linear[k], size_ft, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		  #endif
		 */
  
		// Directly compute the weights; not from precomputed (for small number of nodes, it's very fast)
		generate_conv_weights(conv_weights); 														// calculate the values of the convolution weights (the matrix G_Hat(xi, omega), for xi = (xi_i, xi_j, xi_k), omega = (omega_l, omega_m, omega_n), i,j,k,l,m,n = 0,1,...,N-1) and store the values in conv_weights
		#ifdef FullandLinear																		// only do this FullandLinear was defined
		generate_conv_weights_linear(conv_weights_linear);											// calculate the values of the convolution weights for the linear case (the matrix G_Hat(xi, omega), for xi = (xi_i, xi_j, xi_k), omega = (omega_l, omega_m, omega_n), i,j,k,l,m,n = 0,1,...,N-1) and store the values in conv_weights_linear
		#endif

		MPI_Barrier(MPI_COMM_WORLD);																// set an MPI barrier to ensure that all processes have reached this point before continuing
	}

	char buffer_moment[100], buffer_u[100], buffer_ufull[100], buffer_flags[100],
									buffer_phi[100], buffer_marg[100], buffer_ent[100];				// declare the arrays buffer_moment (to store the name of the file where the moments are printed), buffer_u (to store the name of the file where the solution U is printed), buffer_ufull (to store the name of the file where the solution U is printed in the TwoStream), buffer_flags (to store the flag added to the end of the filenames), buffer_phi (to store the name of the file where the values of phi are printed), buffer_marg (to store the name of the file where the marginals are printed) & buffer_ent (to store the name of the file where the entropy values are printed)

	// EVERY TIME THE CODE IS RUN, CHANGE THE FLAG TO A NAME THAT IDENTIFIES THE CASE RUNNING FOR OR WHAT TIME RUN UP TO:
	sprintf(buffer_flags,"TestFieldFile");																// store the string "OMPfixTest" in buffer_flags
	sprintf(buffer_moment,"Data/Moments_nu%gA%gk%gNx%dLx%gNv%dLv%gSpectralN%ddt%gnT%d_%s.dc",
					nu, A_amp, k_wave, Nx, Lx, Nv, Lv, N, dt, nT, buffer_flags);					// create a .dc file name, located in the directory Data, whose name is Moments_ followed by the values of nu, A_amp, k_wave, Nx, Lx, Nv, Lv, N, dt, nT and the contents of buffer_flags and store it in buffer_moment
	sprintf(buffer_u,"Data/U_nu%gA%gk%gNx%dLx%gNv%dLv%gSpectralN%ddt%gnT%d_%s.dc",
					nu, A_amp, k_wave, Nx, Lx, Nv, Lv, N, dt, nT, buffer_flags);					// create a .dc file name, located in the directory Data, whose name is U_ followed by the values of nu, A_amp, k_wave, Nx, Lx, Nv, Lv, N, dt, nT and the contents of buffer_flags and store it in buffer_u
	sprintf(buffer_ufull,"Data/U2stream_nu%gA%gk%gNx%dLx%gNv%dLv%gSpectralN%ddt%gnT%d_%s.dc",
					nu, A_amp, k_wave, Nx, Lx, Nv, Lv, N, dt, nT, buffer_flags);					// create a .dc file name, located in the directory Data, whose name is U2stream_ followed by the values of nu, A_amp, k_wave, Nx, Lx, Nv, Lv, N, dt, nT and the contents of buffer_flags and store it in buffer_ufull
	sprintf(buffer_marg,"Data/Marginals_nu%gA%gk%gNx%dLx%gNv%dLv%gSpectralN%ddt%gnT%d_%s.dc",
					nu, A_amp, k_wave, Nx, Lx, Nv, Lv, N, dt, nT, buffer_flags);					// create a .dc file name, located in the directory Data, whose name is Marginals_ followed by the values of nu, A_amp, k_wave, Nx, Lx, Nv, Lv, N, dt, nT and the contents of buffer_flags and store it in buffer_moment
	sprintf(buffer_phi,"Data/PhiVals_nu%gA%gk%gNx%dLx%gNv%dLv%gSpectralN%ddt%gnT%d_%s.dc",
					nu, A_amp, k_wave, Nx, Lx, Nv, Lv, N, dt, nT, buffer_flags);					// create a .dc file name, located in the directory Data, whose name is PhiVals_ followed by the values of nu, A_amp, k_wave, Nx, Lx, Nv, Lv, N, dt, nT and the contents of buffer_flags and store it in buffer_moment
	sprintf(buffer_ent,"Data/EntropyVals_nu%gA%gk%gNx%dLx%gNv%dLv%gSpectralN%ddt%gnT%d_%s.dc",
					nu, A_amp, k_wave, Nx, Lx, Nv, Lv, N, dt, nT, buffer_flags);					// create a .dc file name, located in the directory Data, whose name is EntropyVals_ followed by the values of nu, A_amp, k_wave, Nx, Lx, Nv, Lv, N, dt, nT and the contents of buffer_flags and store it in buffer_moment

	#ifdef First																					// only do this if First was defined (setting initial conditions)
		#ifdef Damping																				// only do this if Damping was defined
		SetInit_LD(U);																				// set initial DG solution for Landau Damping. For the first time run t=0, use this to give init solution (otherwise, comment out)
		#endif
		#ifdef TwoStream																			// only do this if TwoStream was defined
		SetInit_LD(U);																				// set initial DG solution for Landau Damping. For the first time run t=0, use this to give init solution (otherwise, comment out)
		#endif
		#ifdef FourHump																				// only do this if FourHump was defined
		SetInit_4H(U);																				// set initial DG solution with the 4Hump IC. For the first time run t=0, use this to give init solution (otherwise, comment out)
		#endif
		#ifdef TwoHump																				// only do this if TwoHump was defined
		SetInit_2H(U);																				// set initial DG solution with the 2Hump IC. For the first time run t=0, use this to give init solution (otherwise, comment out)
		#endif
	#endif

	FILE *fmom, *fu, *fufull, *fmarg, *fphi, *fent;													// declare pointers to the files fmom (which will store the moments), fu (which will store the solution U), fufull (which will store the solution U in the TwoStream case), fmarg (which will store the values of the marginals), fphi (which will store the values of the potential phi) & fent (which will store the values fo the entropy)

	if(myrank_mpi==0)																				// only the process with rank 0 will do this
	{
		#ifdef Damping																				// only do this if Damping was defined
		printf("Damping. %s. nu=%g, A_amp=%g, k_wave=%g, Nx=%d, Lv=%g, Nv=%d, "
				"N=%d, dt=%g, nT=%d\nchunk_Nx=%d, nprocs_Nx=%d\n",
				buffer_flags, nu, A_amp, k_wave, Nx, Lv, Nv, N, dt, nT,chunk_Nx,nprocs_Nx);			// display in the output file that this is the Damping calculation, as well as the contents of the string buffer_flags and the values of nu, A_amp,k_wave, Nx, Lv, Nv, N, dt, nT, chunk_Nx & nprocs_Nx
		#endif
    
		#ifdef TwoStream																			// only do this if TwoStream was defined
		printf("2Stream. %s. with 2Gauss. nu=%g, A_amp=%g, k_wave=%g, Nx=%d, Lv=%g, Nv=%d, "
				"N=%d, dt=%g, nT=%d\n", buffer_flags, nu, A_amp, k_wave, Nx, Lv, Nv, N, dt, nT);	// display in the output file that this is the TwoStream calculation, as well as the contents of the string buffer_flags and the values of nu, A_amp,k_wave, Nx, Lv, Nv, N, dt & nT
		#endif

		#ifdef FourHump																				// only do this if FourHump was defined
		printf("4HumpIC. %s. nu=%g, A_amp=%g, k_wave=%g, Nx=%d, Lv=%g, Nv=%d, "
				"N=%d, dt=%g, nT=%d\nchunk_Nx=%d, nprocs_Nx=%d\n",
				buffer_flags, nu, A_amp, k_wave, Nx, Lv, Nv, N, dt, nT,chunk_Nx,nprocs_Nx);			// display in the output file that this is the calculation with the 4Hump IC, as well as the contents of the string buffer_flags and the values of nu, A_amp,k_wave, Nx, Lv, Nv, N, dt, nT, chunk_Nx & nprocs_Nx
		#endif

		#ifdef TwoHump																				// only do this if TwoHump was defined
		printf("2HumpIC. %s. nu=%g, A_amp=%g, k_wave=%g, Nx=%d, Lv=%g, Nv=%d, "
				"N=%d, dt=%g, nT=%d\nchunk_Nx=%d, nprocs_Nx=%d\n",
				buffer_flags, nu, A_amp, k_wave, Nx, Lv, Nv, N, dt, nT,chunk_Nx,nprocs_Nx);			// display in the output file that this is the calculation with the 2Hump IC, as well as the contents of the string buffer_flags and the values of nu, A_amp,k_wave, Nx, Lv, Nv, N, dt, nT, chunk_Nx & nprocs_Nx
		#endif
	
		#ifdef Second																				// only do this Second was defined (picking up data from a previous run)
		fu=fopen("Data/U_nu0.05A0.2k0.5Nx24Lx12.5664Nv24Lv5.25SpectralN16dt0.01nT1_4HumpMacroTest.dc","r"); 		// set fu to be a file with the name U_nu0.02A0.5k0.785398Nx48Lx8Nv48Lv4.5SpectralN24dt0.004nT500_non_nu002_2.dc, stored in the directory Data, and set the file access mode of fu to r (which means the file must exist already and will be read from)

		fseek(fu ,0 ,SEEK_END);																		// find the final entry of the file fu
	    tp = ftell(fu);																				// use the final entry to determine the size of the file

	    //printf("tp = %d. \n", tp);																// display in the output file the value of tp
	    float NoSol = (float)tp/(size*6*sizeof(double));											// calculate the number of solutions of U that must be in the file
	    printf("no. of set of U = %f \n", NoSol);													// display in the output file the number of solutions of U that must be in the file

		if(tp%(size*6) != 0)																		// check that the number of elements read was a multiple of 6*size
		{
			printf("Error reading file\n");															// if tp was not a multiple of 6*size then display that there was an error
			exit(1);																				// then exit the program
		}
		else
		{
			fseek(fu, tp - size*6*sizeof(double), SEEK_SET);										// find the last solution that was printed in the file
			fread(U,sizeof(double),size*6,fu);														// store the contents of the final solution in the file fu in U, expecting 6*size many entries of datatype double
		}

		fclose(fu);																					// close the file fu
		#endif
      
		fmom=fopen(buffer_moment,"w");																// set fmom to be a file with the name stored in buffer_moment and set the file access mode of fmom to w (which creates an empty file and allows it to be written to)
		fu=fopen(buffer_u, "w");																	// set fu to be a file with the name stored in buffer_u and set the file access mode of fu to w (which creates an empty file and allows it to be written to)
		fmarg=fopen(buffer_marg,"w");																// set fmarg to be a file with the name stored in buffer_marg and set the file access mode of fmarg to w (which creates an empty file and allows it to be written to)
		fphi=fopen(buffer_phi,"w");																	// set fphi to be a file with the name stored in buffer_phi and set the file access mode of fphi to w (which creates an empty file and allows it to be written to)
		fent=fopen(buffer_ent,"w");																	// set fent to be a file with the name stored in buffer_ent and set the file access mode of fent to w (which creates an empty file and allows it to be written to)

		FindNegVals(U, fNegVals, fAvgVals);															// find out in which cells the approximate solution goes negative and record it in fNegVals

		/* DEBUG TEST
		for(int i=0;i<Nx;i++)
		{
			printf("i = %d: ", i);
			for(int j1=0;j1<Nv;j1++)
			{
				for(int j2=0;j2<Nv;j2++)
				{
					for(int j3=0;j3<Nv;j3++)
					{
						printf("%g ", fEquiVals[5*(5*(5*(5*(i*size_v + j1*Nv*Nv+j2*Nv+j3)+2) + 2) + 2) + 2]);
					}
				}
			}
			printf("\n");
		}
		*/

		mass=computeMass(U);																		// set mass to the value calculated through computeMass, for the solution f(x,v,t) at the current time t, using its DG coefficients stored U
		computeMomentum(U, a);																		// calculate the momentum for the solution f(x,v,t) at the current time t, using its DG coefficients stored U, and store it in a
		KiE=computeKiE(U);																			// set KiE to the value calculated through computeKiE, for the solution f(x,v,t) at the current time t, using its DG coefficients stored U
		EleE=computeEleE(U);																		// set EleE to the value calculated through computeEleE, for the solution f(x,v,t) at the current time t, using its DG coefficients stored U
		tmp = sqrt(EleE);																			// set tmp to the square root of EleE
		ent1 = computeEntropy(U);																	// set ent1 the value calculated through computeEntropy
		l_ent1 = log(fabs(ent1));																	// set l_ent1 to the log of ent1
		ll_ent1 = log(fabs(l_ent1));																// set ll_ent1 to the log of l_ent1
		printf("step #0: %11.8g  %11.8g  %11.8g  %11.8g  %11.8g  %11.8g  %11.8g  %11.8g %11.8g %11.8g \n",
				mass, a[0], a[1], a[2], KiE, EleE, tmp, log(tmp), KiE+EleE, ent1);					// display in the output file that this is step 0 (so these are the initial conditions), then the mass, 3 components of momentum, kinetic energy, electric energy, sqrt(electric energy), log(sqrt(electric energy)), total energy & entropy
		fprintf(fmom, "%11.8g %11.8g %11.8g  %11.8g  %11.8g  %11.8g  %11.8g  %11.8g  %11.8g \n",
				mass, a[0], a[1], a[2], KiE, EleE, tmp, log(tmp), KiE+EleE);						// in the file tagged as fmom, print the initial mass, 3 components of momentum, kinetic energy, electric energy, sqrt(electric energy), log(sqrt(electric energy)) & total energy
		fprintf(fent, "%11.8g %11.8g %11.8g \n", ent1, l_ent1, ll_ent1);							// in the file tagged as fent, print the entropy, its log and the log of that

		KiEratio = computeKiEratio(U, fNegVals);													// compute the ratio of the kinetic energy where f is negative to that where it is positive and store it in KiEratio
		printf("Kinetic Energy Ratio = %g\n", KiEratio);											// print the ratio of the kinetic energy where f is negative to that where it is positive

		//fufull=fopen("Data/U_nu0.02A0.5k1.5708Nx48Lx4Nv32Lv4SpectralN24dt0.004_non_nu002_time15s.dc", "w");
		//fprintf(fmom, "%11.8g  %11.8g\n", EleE, log(tmp));
		/*#ifdef TwoStream
		fufull=fopen(buffer_ufull, "w");
		for(k=0;k<size;k++){
		for(l=0;l<6;l++)fprintf(fufull, "%g ", U[k*6+l]);
		}
		fprintf(fufull,"\n\n");
		#endif*/

		PrintMarginalLoc(fmarg);																	// print the values of x & v1 that the marginal will be evaluated at in the file tagged as fmarg
		PrintMarginal(U, fmarg);																	// print the marginal distribution for the initial condition, using the DG coefficients in U, in the file tagged as fmarg
	}
  
	MPI_Bcast(U, size*6, MPI_DOUBLE, 0, MPI_COMM_WORLD);   											// send the contents of U, which will be 6*size entries of datatype MPI_DOUBLE, from the process with rank 0 to all processes, using the communicator MPI_COMM_WORLD
	MPI_Barrier(MPI_COMM_WORLD);																	// set an MPI barrier to ensure that all processes have reached this point before continuing
  
	MPIt1 = MPI_Wtime();																			// set MPIt1 to the current time in the MPI process
	while(t < nT) 																					// if t < nT (i.e. not yet reached the final timestep), perform time-splitting to first advect the particle through the collisionless step and then perform one space homogeneous collisional step)
	{
		RK3(U); 																					// Use RK3 to perform one timestep of the collisionless problem

		if(nu > 0.)
		{
			setInit_spectral(U, f); 																// Take the coefficient of the DG solution from the advection step, and project them onto the grid used for the spectral method to perform the collision step

			for(l=chunk_Nx*myrank_mpi;l<chunk_Nx*(myrank_mpi+1) && l<Nx;l++) 						// divide the number of discretised space points equally over the number of MPI processes, so that each process receives a different chunk of space to work on
			{
				#ifdef FullandLinear																// only do this if FullandLinear was defined
				ComputeQ(f[l%chunk_Nx], qHat, conv_weights, qHat_linear, conv_weights_linear);		// using the coefficients of the current solution stored in f (but only for the chunk of space being taken care of by the current MPI process), calculate the Fourier tranform of Q(f,f) using conv_weights for the weights in the convolution in the full part of Q & conv_weights_linear in the convolution in the linear part of Q, then store the results of each Fourier transform in qHat & qHat_linear, respectively
				conserveAllMoments(qHat, qHat_linear);												// perform the explicit conservation calculation
				RK4(f[l%chunk_Nx], l, qHat, conv_weights, qHat_linear, conv_weights_linear,
						U, Utmp_coll);																// advance to the next time step in the collisional problem using RK4 at the given space-step l, taking the current solution stored in f (but only for the chunk of space being taken care of by the current MPI process), as well as qHat, conv_weights, qHat_linear & conv_weights_linear (to allow more Fourier transforms of Q to be made), storing the output partially in U and partially in Utmp_coll
				#else																				// otherwise, if FullandLinear was not defined...
				ComputeQ(f[l%chunk_Nx], qHat, conv_weights);										// using the coefficients of the current solution stored in f (but only for the chunk of space being taken care of by the current MPI process), calculate the Fourier tranform of Q(f,f) using conv_weights for the weights in the convolution, then store the results of the Fourier transform in qHat
				conserveAllMoments(qHat);															// perform the explicit conservation calculation
				RK4(f[l%chunk_Nx], l, qHat, conv_weights, U, Utmp_coll);							// advance to the next time step in the collisional problem using RK4 at the given space-step l, taking the current solution stored in f (but only for the chunk of space being taken care of by the current MPI process), as well as qHat & conv_weights (to allow more Fourier transforms of Q to be made), storing the output partially in U and partially in Utmp_coll
				#endif
			}

			MPI_Barrier(MPI_COMM_WORLD);															// set an MPI barrier to ensure that all processes have reached this point before continuing
			if(myrank_mpi == 0) 																	// only the process with rank 0 will do this
			{
				// TRANSFER CONTENTS OF THE dU (Utmp_coll) THAT HAVE BEEN COMPUTED INTO U1 (U):
				for(l=0;l<chunk_Nx;l++) 															// cycle through all space-steps stored in the first chunk of the space interval (which is stored on the process with rank 0)
				{
					for(k=0;k<size_v;k++)															// cycle through all size_v (= Nv^3) many velocity-steps (which will exist for each space-step)
					{
						k_v = l*size_v + k;															// set k_v to be the value associated with the k-th velocity-step for the l-th space-step
						U[k_v*6+0] = Utmp_coll[k_v*5];												// set the 6*k_v-th entry of U to the 5*k_v-th entry of Utmp_coll
						U[k_v*6+5] = Utmp_coll[k_v*5+4]; 											// set the (6*k_v + 5)-th entry of U to the (5*k_v + 4)-th entry of Utmp_coll
						U[k_v*6+2] = Utmp_coll[k_v*5+1];  											// set the (6*k_v + 2)-th entry of U to the (5*k_v + 1)-th entry of Utmp_coll
						U[k_v*6+3] = Utmp_coll[k_v*5+2];	 										// set the (6*k_v + 3)-th entry of U to the (5*k_v + 2)-th entry of Utmp_coll
						U[k_v*6+4] = Utmp_coll[k_v*5+3];											// set the (6*k_v + 4)-th entry of U to the (5*k_v + 3)-th entry of Utmp_coll
					}
				}
				// RECEIVE FROM ALL OTHER PROCESSES CONSECUTIVELY TO ENSURE THE WEIGHTS ARE STORED IN THE FILE U CONSECUTIVELY:
				for(i=1;i<nprocs_Nx;i++)															// store the DG coefficients of the current solution in U that were calculated by the remaining processes (with ranks i = 1, 2, ..., nprocs_Nx-1) for their corresponding chunk of space
				{
					MPI_Recv(output_buffer, chunk_Nx*size_v*5, MPI_DOUBLE, i, i,
							MPI_COMM_WORLD, &status);											 	// receive a message of 5*chunk_Nx_size_v entries of datatype MPI_DOUBLE from the process with rank i (storing the i-th space chunk of U, containing the DG coefficients calculate on the processor with corresponding rank), storing the data in output_buffer, with tag i in the communicator MPI_COMM_WORLD, storing the status of the receive in status
					for(l=0;l<chunk_Nx;l++)															// cycle through all space-steps stored in the chunk of the space interval dealt with by the i-th process
					{
						if((chunk_Nx*i+l)<Nx)														// ensure that the space-step being dealt with exists
						{
							for(k=0;k<size_v;k++)													// cycle through all size_v (= Nv^3) many velocity-steps (which will exist for each space-step)
							{
								k_v = (chunk_Nx*i+l)*size_v + k; 									// set k_v to be the value associated with the k-th velocity-step for the (chunk_Nx*i + l)-th space-step (which is the l-th space-step in the current chunk)
								k_local = l*size_v + k;												// set k_local to be the value associated with the k-th veolcity-step for the l-th local space step (which is the current one, as indicated in the previous line)

								// STORE CONTENTS OF THE RECEIVE BUFFER IN THE CORRECT PORTION OF U TO ADD THIS PART OF THE SOLUTION:
								U[k_v*6+0] = output_buffer[k_local*5];								// set the 6*k_v-th entry of U to the 5*k_local-th entry of Utmp_coll
								U[k_v*6+5] = output_buffer[k_local*5+4];							// set the (6*k_v + 5)-th entry of U to the (5*k_local + 4)-th entry of Utmp_coll
								U[k_v*6+2] = output_buffer[k_local*5+1];  							// set the (6*k_v + 2)-th entry of U to the (5*k_local + 1)-th entry of Utmp_coll
								U[k_v*6+3] = output_buffer[k_local*5+2];	 						// set the (6*k_v + 3)-th entry of U to the (5*k_local + 2)-th entry of Utmp_coll
								U[k_v*6+4] = output_buffer[k_local*5+3];							// set the (6*k_v + 4)-th entry of U to the (5*k_local + 3)-th entry of Utmp_coll
							}
						}
					}
				}
			}
			else 																					// the remaining processes, with rank 1, 2, ..., nprocs_Nx-1 will do this
			{
				if(myrank_mpi<nprocs_Nx)
				{
					MPI_Send(Utmp_coll, chunk_Nx*size_v*5, MPI_DOUBLE, 0, myrank_mpi,
							MPI_COMM_WORLD);														// send the contents of Utmp_coll, which will be 5*chunk_Nx*size_v entries of datatype MPI_DOUBLE, to the process with rank 0, tagged with the rank of the current process, via the MPI_COMM_WORLD communicator
				}
			}
		    MPI_Bcast(U, size*6, MPI_DOUBLE, 0, MPI_COMM_WORLD);    								// send the contents of U, from the process with rank 0, which contains 6*size entries of datatype MPI_DOUBLE, to all processes via the communicator MPI_COMM_WORLD (so that all processes have the coefficients of the DG approximation to f at the current time-step for the start of the next calculation)
		}
   
		if(myrank_mpi==0)																			// only the process with rank 0 will do this
		{
			FindNegVals(U, fNegVals, fAvgVals);																// find out in which cells the approximate solution goes negative and record it in fNegVals

			mass=computeMass(U);																	// set mass to the value calculated through computeMass, for the solution f(x,v,t) at the current time t, using its DG coefficients stored U
			computeMomentum(U, a);																	// calculate the momentum for the solution f(x,v,t) at the current time t, using its DG coefficients stored U, and store it in a
			KiE=computeKiE(U);																		// set KiE to the value calculated through computeKiE, for the solution f(x,v,t) at the current time t, using its DG coefficients stored U
			EleE=computeEleE(U);																	// set EleE to the value calculated through computeEleE, for the solution f(x,v,t) at the current time t, using its DG coefficients stored U
			tmp = sqrt(EleE);																		// set tmp to the square root of EleE
			ent1 = computeEntropy(U);																// set ent1 the value calculated through computeEntropy
			l_ent1 = log(fabs(ent1));																// set l_ent1 to the log of ent1
			ll_ent1 = log(fabs(l_ent1));															// set ll_ent1 to the log of l_ent1
			printf("step %d: %11.8g  %11.8g  %11.8g  %11.8g  %11.8g  %11.8g  %11.8g  %11.8g %11.8g %11.8g \n",
					t+1, mass, a[0], a[1], a[2], KiE, EleE, tmp, log(tmp), KiE+EleE, ent1);			// display in the output file that this is step t+1, then the mass, 3 components of momentum, kinetic energy, electric energy, sqrt(electric energy), log(sqrt(electric energy)), total energy & entropy
			fprintf(fmom, "%11.8g %11.8g %11.8g  %11.8g  %11.8g  %11.8g  %11.8g  %11.8g  %11.8g \n",
					mass, a[0], a[1], a[2], KiE, EleE, tmp, log(tmp), KiE+EleE);					// in the file tagged as fmom, print the initial mass, 3 components of momentum, kinetic energy, electric energy, sqrt(electric energy), log(sqrt(electric energy)) & total energy
			fprintf(fent, "%11.8g %11.8g %11.8g \n", ent1, l_ent1, ll_ent1);						// in the file tagged as fent, print the entropy, its log and the log of that

			KiEratio = computeKiEratio(U, fNegVals);												// compute the ratio of the kinetic energy where f is negative to that where it is positive and store it in KiEratio
			printf("Kinetic Energy Ratio = %g\n", KiEratio);										// print the ratio of the kinetic energy where f is negative to that where it is positive

			//fprintf(fmom, "%11.8g  %11.8g\n", EleE, log(tmp));
			/*#ifdef TwoStream
			  if(t==20 || t==50 || t==80 || t==250)
			  {
				  for(k=0;k<size;k++)
				  {
						for(l=0;l<6;l++)
						{
							fprintf(fufull, "%g ", U[k*6+l]);
						}
				  }
				  fprintf(fufull,"\n\n");
			  }
			  #endif*/
      

	    	//if(t%400==0)fwrite(U,sizeof(double),size*6,fu);
			if(t%20==0)
			{
				PrintMarginal(U, fmarg);															// print the marginal distribution for the initial condition, using the DG coefficients in U, in the file tagged as fmarg
			}
		}
	
		t++;																						// increment t by one
	
		MPI_Barrier(MPI_COMM_WORLD);																// set an MPI barrier to ensure that all processes have reached this point before continuing
	}
  
	MPIelapsed = MPI_Wtime() - MPIt1;																// set MPIelapsed to the current time minus MPIt1 to calculate how long nT time-steps took
	if(myrank_mpi==0)																				// only the process with rank 0 will do this
	{
		printf("time duration for %d time steps is %gs\n",nT, MPIelapsed);							// display in the output file how long it took to calculate nT time-steps
    
		fwrite(U,sizeof(double),size*6,fu);															// write the coefficients of the DG approximation at the end, stored in U, which is 6*size entires, each of the size of a double datatype, in the file tagged as fu
		PrintPhiVals(U, fphi);																		// print the values of the potential in the file tagged as filephi at the given timestep
	
		fclose(fu);  																				// remove the tag fu to close the file
	}
	MPI_Barrier(MPI_COMM_WORLD);																	// set an MPI barrier to ensure that all processes have reached this point before continuing
  
	if(myrank_mpi==0)																				// only the process with rank 0 will do this
	{
		fclose(fmom);  																				// remove the tag fmom to close the file
		fclose(fmarg);  																			// remove the tag fmarg to close the file
		fclose(fphi);  																				// remove the tag fphi to close the file
		fclose(fent);  																				// remove the tag fent to close the file
		#ifdef TwoStream																			// only do this if TwoStream was defined
		//fclose(fufull);
		#endif
	}
	if(nu > 0.)
	{
		free(C1); free(C2); free(v); free(eta); free(wtN);											// delete the dynamic memory allocated for C1, C2, v, eta & wtN
		free(f); free(conv_weights); free(output_buffer); 											// delete the dynamic memory allocated for f, conv_weights, output_buffer
		fftw_free(temp); fftw_free(qHat);															// delete the dynamic memory allocated for temp & qhat
		fftw_free(Q1_fft); fftw_free(Q2_fft); fftw_free(Q3_fft); fftw_free(fftOut); fftw_free(fftIn); // delete the dynamic memory allocated for Q1_fft, Q2_fft, Q3_fft, fftOut & fftIn
		free(Q);free(f1);free(Q1); free(Utmp_coll);// free(f2); free(f3);//free(Q3);				// delete the dynamic memory allocated for Q, f1, Q1 & Utmp_coll
		#ifdef FullandLinear																		// only do this if FullandLinear is defined
		fftw_free(qHat_linear); fftw_free(Q1_fft_linear); 											// delete the dynamic memory allocated for qHat_linear & Q1_fft_linear
		fftw_free(Q2_fft_linear); fftw_free(Q3_fft_linear); 										// delete the dynamic memory allocated for Q2_fft_linear & Q3_fft_linear
		free(conv_weights_linear);																	// delete the dynamic memory allocated for conv_weights_linear
		#endif
	}
	free(output_buffer_vp);																			// delete the dynamic memory allocated for output_buffer_vp
	free(U); free(U1); free(Utmp); // free(H);														// delete the dynamic memory allocated for U, U1 & Utmp
	free(cp); free(intE); free(intE1); free(intE2);													// delete the dynamic memory allocated for cp, intE, intE1 & inteE2

	free(fNegVals); free(fAvgVals);	free (fEquiVals);												// delete the dynamic memory allocated for fNegVals, fAvgVals & fEquiVals
  
	iparse.Close();																					// close the input file

	MPI_Finalize();																					// ensure that MPI exits cleanly
	return 0;																						// return 0, since main is of type int (and this shows the program completed correctly)
}
