/* This is the source file which contains the subroutines necessary for solving the space homogeneous,
 * collision problem resulting from time-splitting, including FFT routines.
 *
 * Functions included: S1hat, S233hat, S213hat, gHat3, gHat3_linear, generate_conv_weights,
 * generate_conv_weights_linear, fft3D, ifft3D, FS, ComputeQ, IntModes, ProjectedNodeValue, RK4
 *
 */

#include "collisionRoutines_1.h"																		// collisionRoutines_1.h is where the prototypes for the functions contained in this file are declared

extern fftw_plan p_forward; 
extern fftw_plan p_backward; 
extern fftw_complex *temp;

double IntM[10];																					// declare an array IntM to hold 10 double variables
#pragma omp threadprivate(IntM)																		// start the OpenMP parallel construct to start the threads which will run in parallel, passing IntM to each thread as private variables which will have their contents deleted when the threads finish (doesn't seem to be doing anything since no {} afterwards???)

double S1hat(double ki1,double ki2,double ki3)
{
  if(ki1==0. && ki2==0. && ki3==0.) return sqrt(1./(2*PI))*R_v*R_v;
  else return sqrt(2.0/PI)*(1-cos(R_v*sqrt(ki1*ki1+ki2*ki2+ki3*ki3)))/(ki1*ki1+ki2*ki2+ki3*ki3);
}

double S233hat(double ki1, double ki2, double ki3)
{
  double r=sqrt(ki1*ki1+ki2*ki2+ki3*ki3);
  if(r==0.) return sqrt(1./(2.*PI))*R_v*R_v/3.;
  else return  sqrt(2./PI)*((ki1*ki1+ki2*ki2)*(R_v*r-sin(R_v*r))/(R_v*r)-ki3*ki3*(R_v*r+R_v*r*cos(R_v*r)-2.*sin(R_v*r))/(R_v*r))/pow(r,4.);
}

double S213hat(double ki1, double ki2, double ki3)
{
  double r=sqrt(ki1*ki1+ki2*ki2+ki3*ki3);
  if(ki1==0. || ki3==0.) return 0.;
  else return -sqrt(2/PI)*ki1*ki3*(2.*R_v*r+R_v*r*cos(R_v*r)-3.*sin(R_v*r))/(R_v*pow(r,5.));
}

double gHat3(double eta1, double eta2, double eta3, double ki1, double ki2, double ki3 ) 
{
	double result = 0.;
	double ki[3]={ki1,ki2,ki3}, zeta[3]={eta1,eta2,eta3}; // ki=w, zeta=xi in the notes
	double Shat[3][3];
	double r=sqrt(ki1*ki1+ki2*ki2+ki3*ki3);
	//double darg[3][2];
	int i,j;
	
	/*Shat[0][0]=S1hat(-ki1,-ki2,-ki3)-S233hat(-ki2,-ki3,-ki1)-S1hat(eta1-ki1,eta2-ki2,eta3-ki3)+S233hat(eta2-ki2,eta3-ki3,eta1-ki1);
	Shat[1][1]=S1hat(-ki1,-ki2,-ki3)-S233hat(-ki1,-ki3,-ki2)-S1hat(eta1-ki1,eta2-ki2,eta3-ki3)+S233hat(eta1-ki1,eta3-ki3,eta2-ki2);
	Shat[2][2]=S1hat(-ki1,-ki2,-ki3)-S233hat(-ki1,-ki2,-ki3)-S1hat(eta1-ki1,eta2-ki2,eta3-ki3)+S233hat(eta1-ki1,eta2-ki2,eta3-ki3);
	Shat[0][1]=-S212hat(-ki1,-ki2,-ki3)+S212hat(eta1-ki1,eta2-ki2,eta3-ki3);
	Shat[0][2]=-S212hat(-ki1,-ki3,-ki2)+S212hat(eta1-ki1,eta3-ki3,eta2-ki2);
	Shat[1][2]=-S212hat(-ki2,-ki3,-ki1)+S212hat(eta2-ki2,eta3-ki3,eta1-ki1);
	Shat[1][0]=Shat[0][1]; Shat[2][0]=Shat[0][2]; Shat[2][1]=Shat[1][2]; */

	Shat[0][0]=S1hat(ki1,ki2,ki3)-S233hat(ki2,ki3,ki1);
	Shat[1][1]=S1hat(ki1,ki2,ki3)-S233hat(ki1,ki3,ki2);
	Shat[2][2]=S1hat(ki1,ki2,ki3)-S233hat(ki1,ki2,ki3);
	Shat[0][1]=-S213hat(ki1,ki3,ki2);
	Shat[0][2]=-S213hat(ki1,ki2,ki3);
	Shat[1][2]=-S213hat(ki2,ki1,ki3);
	Shat[1][0]=Shat[0][1]; Shat[2][0]=Shat[0][2]; Shat[2][1]=Shat[1][2];
	
	/*darg[0][0]=ki[0];darg[0][1]=sqrt(ki[1]*ki[1]+ki[2]*ki[2]);
	darg[1][0]=ki[1];darg[1][1]=sqrt(ki[2]*ki[2]+ki[0]*ki[0]);
	darg[2][0]=ki[2];darg[2][1]=sqrt(ki[0]*ki[0]+ki[1]*ki[1]);
	for(i=0;i<3;i++)result += 4*gauss_legendre(GL, Divg0, darg, 0, 1)*zeta[i];*/

	for(i=0;i<3;i++){
	  for(j=0;j<3;j++){
	    //result += Shat[i][j]*(2.*ki[j]-zeta[j])*zeta[i];
	    result += Shat[i][j]*(zeta[i]-ki[i])*(zeta[j]-ki[j]);
	  }
	}	
	if(r==0.)result = -result;
	else result = (sqrt(8./PI))*(R_v*r-sin(R_v*r))/(R_v*r) - result;
  	return result;	
}

double gHat3_linear(double eta1, double eta2, double eta3, double ki1, double ki2, double ki3 ) 
{
	double result = 0.;
	double ki[3]={ki1,ki2,ki3}, zeta[3]={eta1,eta2,eta3}; // ki=w, zeta=xi in the notes
	double Shat[3][3];
	double r=sqrt(ki1*ki1+ki2*ki2+ki3*ki3);
	//double darg[3][2];
	int i,j;

	Shat[0][0]=S1hat(ki1,ki2,ki3)-S233hat(ki2,ki3,ki1);
	Shat[1][1]=S1hat(ki1,ki2,ki3)-S233hat(ki1,ki3,ki2);
	Shat[2][2]=S1hat(ki1,ki2,ki3)-S233hat(ki1,ki2,ki3);
	Shat[0][1]=-S213hat(ki1,ki3,ki2);
	Shat[0][2]=-S213hat(ki1,ki2,ki3);
	Shat[1][2]=-S213hat(ki2,ki1,ki3);
	Shat[1][0]=Shat[0][1]; Shat[2][0]=Shat[0][2]; Shat[2][1]=Shat[1][2];

	for(i=0;i<3;i++){
	  for(j=0;j<3;j++){
	    //result += Shat[i][j]*(2.*ki[j]-zeta[j])*zeta[i];
	    result += Shat[i][j]*zeta[i]*(zeta[j]-ki[j]);
	  }
	}	
	result = -result;
  	return result;	
}

/*#ifdef UseMPI
void generate_conv_weights(double **conv_weights)
{
  int t, i, j, k, l, m, n;
  
   #pragma omp parallel for  private(i,j,k,l,m,n) shared(conv_weights)
 for(t=chunksize_ft*myrank_mpi;t<chunksize_ft*(myrank_mpi+1);t++){   
      k = t % N;
      j = ((t-k)/N) % N;
      i = (t - k - N*j)/(N*N);
      for(l=0;l<N;l++){
	for(m=0;m<N;m++){
	  for(n=0;n<N;n++) {
	    conv_weights[t%chunksize_ft][n + N*(m + N*l)] =gHat3(eta[i], eta[j], eta[k], eta[l], eta[m], eta[n]); 
	  }
	}
      }      
  }
}
#else */
void generate_conv_weights(double **conv_weights)
{
  int i, j, k, l, m, n;
   #pragma omp parallel for private(i,j,k,l,m,n) shared(conv_weights)
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      for(k=0;k<N;k++){ 
	for(l=0;l<N;l++){
	  for(m=0;m<N;m++){
	    for(n=0;n<N;n++) {
	     conv_weights[k + N*(j + N*i)][n + N*(m + N*l)] = gHat3(eta[i], eta[j], eta[k], eta[l], eta[m], eta[n]); // in the notes, correspondingly, (i,j,k)-kxi, (l,m,n)-w	      	      
	    }
	  }
	}
      }
    }
  }
}

void generate_conv_weights_linear(double **conv_weights_linear)
{
  int i, j, k, l, m, n;
   #pragma omp parallel for private(i,j,k,l,m,n) shared(conv_weights_linear)
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      for(k=0;k<N;k++){ 
	for(l=0;l<N;l++){
	  for(m=0;m<N;m++){
	    for(n=0;n<N;n++) {
	     conv_weights_linear[k + N*(j + N*i)][n + N*(m + N*l)] = gHat3_linear(eta[i], eta[j], eta[k], eta[l], eta[m], eta[n]); // in the notes, correspondingly, (i,j,k)-kxi, (l,m,n)-w	      	      
	    }
	  }
	}
      }
    }
  }
}
//#endif
/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

/*
function fft3D
--------------
Computes the fourier transform of in, and adjusts the coefficients based on our v, eta
*/
void fft3D(fftw_complex *in, fftw_complex *out)
{
  int i, j, k, index;
  double sum;
  
  //shift the 'v' terms in the exponential to reflect our velocity domain
  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      for(k=0;k<N;k++)
	{

	  index = k + N*(j + N*i);
	  sum = ((double)i + (double)j + (double)k)*L_eta*h_v;

	  //h_v correspond to the velocity space scaling - ensures that the FFT is properly scaled since fftw does no scaling at all
	  temp[index][0] = scale3*h_v*h_v*h_v*wtN[i]*wtN[j]*wtN[k]*(cos(sum)*in[index][0] - sin(sum)*in[index][1]);
	  temp[index][1] = scale3*h_v*h_v*h_v*wtN[i]*wtN[j]*wtN[k]*(cos(sum)*in[index][1] + sin(sum)*in[index][0]);
	}
  //computes fft
  fftw_execute(p_forward);
  //fftwnd_threads_one(nThreads, p_forward, temp, NULL); /*FFTW Library*/

  //shifts the 'eta' terms to reflect our fourier domain
  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      for(k=0;k<N;k++)
	{
	  index = k + N*(j + N*i);
	  sum = L_v*(eta[i] + eta[j] + eta[k]);
	  
	  out[index][0] = ( cos(sum)*temp[index][0] - sin(sum)*temp[index][1]);
	  out[index][1] = ( cos(sum)*temp[index][1] + sin(sum)*temp[index][0]);
	}
  
}

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
/*
function fft3D
--------------
Computes the inverse fourier transform of in, and adjusts the coefficients based on our v, eta
*/

void ifft3D(fftw_complex *in, fftw_complex *out)
{
  int i, j, k, index;
  double sum, numScale = scale3;//= pow((double)N, -3.0);

  //shifts the 'eta' terms to reflect our fourier domain
  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      for(k=0;k<N;k++)
	{
	  index = k + N*(j + N*i);
	  sum = -( ( (double)i + (double)j + (double)k )*L_v*h_eta );
	  
	  //h_eta ensures FFT is scaled correctly, since fftw does no scaling at all
	  temp[index][0] = h_eta*h_eta*h_eta*wtN[i]*wtN[j]*wtN[k]*(cos(sum)*in[index][0] - sin(sum)*in[index][1]);
	  temp[index][1] = h_eta*h_eta*h_eta*wtN[i]*wtN[j]*wtN[k]*(cos(sum)*in[index][1] + sin(sum)*in[index][0]);
	}
  //compute IFFT
  fftw_execute(p_backward);
  //fftwnd_threads_one(nThreads, p_backward, temp, NULL); /*FFTW Library*/

  //shifts the 'v' terms to reflect our velocity domain
  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      for(k=0;k<N;k++)
	{
	  index = k + N*(j + N*i);
	  sum = -( L_eta*(v[i] + v[j] + v[k])  );
	  
	  out[index][0] = (cos(sum)*temp[index][0] - sin(sum)*temp[index][1])*numScale;
	  out[index][1] = (cos(sum)*temp[index][1] + sin(sum)*temp[index][0])*numScale;
	}
  
}

void FS(fftw_complex *in, fftw_complex *out) // compute the Fourier series approximation of f (out), through fhat (in)
{
  int i, j, k, index;
  double sum;//= pow((double)N, -3.0);

  //shifts the 'eta' terms to reflect our fourier domain
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      for(k=0;k<N;k++){
	  index = k + N*(j + N*i);
	  sum = -( ( (double)i + (double)j + (double)k )*L_v*h_eta );
	  
	  //h_eta ensures FFT is scaled correctly, since fftw does no scaling at all
	  temp[index][0] = (cos(sum)*in[index][0] - sin(sum)*in[index][1]);
	  temp[index][1] = (cos(sum)*in[index][1] + sin(sum)*in[index][0]);
      }
    }
  }
  //compute IFFT
  fftw_execute(p_backward);
  //fftwnd_threads_one(nThreads, p_backward, temp, NULL); /*FFTW Library*/

  //shifts the 'v' terms to reflect our velocity domain
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      for(k=0;k<N;k++){
	  index = k + N*(j + N*i);
	  sum = -( L_eta*(v[i] + v[j] + v[k])  );
	  
	  out[index][0] = (cos(sum)*temp[index][0] - sin(sum)*temp[index][1])/scaleL/scale3;
	  out[index][1] = (cos(sum)*temp[index][1] + sin(sum)*temp[index][0])/scaleL/scale3;
      }
    }
  }
  
}

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

/*
function ComputeQ
-----------------
The main function for calculating the collision effects
 */

#ifdef MPI_parallelcollision
void ComputeQ(double *f, fftw_complex *qHat, double **conv_weights)
{
  int t, i, j, k, l, m, n, x, y, z;
  int start_i, start_j, start_k, end_i, end_j, end_k;
  //fftw_complex *fftIn, *fftOut;
  double tempD, tmp0, tmp1;
  double prefactor = h_eta*h_eta*h_eta; // we dont have scale3 here.

  //fftIn = (fftw_complex *)fftw_malloc(N*N*N*sizeof(fftw_complex));
  //fftOut = (fftw_complex *)fftw_malloc(N*N*N*sizeof(fftw_complex));  

  for(i=0;i<size_ft;i++){
    fftIn[i][0] = f[i];
    fftIn[i][1] = 0.;
  }
  fft3D(fftIn, fftOut);
  
  #pragma omp parallel for schedule(dynamic) private(j,k,l,m,n,x,y,z,start_i,start_j,start_k,end_i,end_j,end_k,tempD, tmp0, tmp1) shared(qHat, conv_weights)
  for(t=chunksize_ft*myrank_mpi;t<chunksize_ft*(myrank_mpi+1);t++){   
      t_local = t%chunksize_ft;
      k = t % N;
      j = ((t-k)/N) % N;
      i = (t - k - N*j)/(N*N);
	  //figure out the windows for the convolutions (i.e. where xi(l) and eta(i)-xi(l) are in the domain)
	  if( i < N/2 ) {
	    start_i = 0;
	    end_i = i + N/2 + 1; 
	  }
	  else {
	    start_i = i - N/2 + 1;
	    end_i = N;
	  }
	  
	  if( j < N/2 ) {
	    start_j = 0;
	    end_j = j + N/2 + 1; 
	  }
	  else {
	    start_j = j - N/2 + 1;
	    end_j = N;
	  }
	  
	  if( k < N/2 ) {
	    start_k = 0;
	    end_k = k + N/2 + 1; 
	  }
	  else {
	    start_k = k - N/2 + 1;
	    end_k = N;
	  }
	  tmp0=0.; tmp1=0.;
	  for(l=start_i;l<end_i;l++) {
	    for(m=start_j;m<end_j;m++) {
	      for(n=start_k;n<end_k;n++)
		{
		  x = i + N/2 - l;
		  y = j + N/2 - m;
		  z = k + N/2 - n;
		  //get convolution weight
		  tempD = conv_weights[t_local][n + N*(m + N*l)];                 
		  tmp0 += prefactor*wtN[l]*wtN[m]*wtN[n]*tempD*(fftOut[n + N*(m + N*l)][0]*fftOut[z + N*(y + N*x)][0] - fftOut[n + N*(m + N*l)][1]*fftOut[z + N*(y + N*x)][1]);
		  
		  tmp1 += prefactor*wtN[l]*wtN[m]*wtN[n]*tempD*(fftOut[n + N*(m + N*l)][0]*fftOut[z + N*(y + N*x)][1] + fftOut[n + N*(m + N*l)][1]*fftOut[z + N*(y + N*x)][0]);
		}
	    }
	  }
	 qHat[t_local][0] = tmp0;  
	 qHat[t_local][1] = tmp1;
    }  

 // fftw_free(fftIn);
 // fftw_free(fftOut);

}
#endif


void IntModes(int k1, int k2,  int k3, int j1, int j2, int j3, double *result) // \int_Ij exp(i\xi_k \cdot v) \phi_l(v) dv; l=0..4: 1, (v1-w_j1)/dv, (v2-w_j2)/dv, (v3-w_j3)/dv, square sum of last three components
{
  double tmp_re, tmp_im, tmp1_re, tmp1_im, tmp2_re, tmp2_im, tmp3_re, tmp3_im, tem1_re, tem1_im, tem2_re, tem2_im, tem3_re, tem3_im, tp1_re, tp1_im, tp2_re, tp2_im, tp3_re, tp3_im, a_re, a_im;

  double v1, v2, v3;																	// define v1, v2 & v3 (where v_k is the center of the current velocity cell in the kth direction)
  double v1_l, v1_r, v2_l, v2_r, v3_l, v3_r;											// define v1_l, v1_r, v2_l, v2_r, v3_l & v3_r (where vk_l & vk_r are the values of v_(j_k) at the left and right edges of the velocity cell j in the kth direction)
  v1 = Gridv((double)j1);																		// calculate v1 by the formula -Lv+(j1+0.5)*dv
  v2 = Gridv((double)j2);																		// calculate v2 by the formula -Lv+(j2+0.5)*dv
  v3 = Gridv((double)j3);																		// calculate v3 by the formula -Lv+(j3+0.5)*dv
  v1_l = Gridv(j1-0.5);																	// calculate v1_l by the formula -Lv+j1*dv
  v1_r = Gridv(j1+0.5);																	// calculate v1_r by the formula -Lv+(j1+1)*dv
  v2_l = Gridv(j2-0.5);																	// calculate v2_l by the formula -Lv+j2*dv
  v2_r = Gridv(j2+0.5);																	// calculate v2_r by the formula -Lv+(j2+1)*dv
  v3_l = Gridv(j3-0.5);																	// calculate v3_l by the formula -Lv+j3*dv
  v3_r = Gridv(j3+0.5);																	// calculate v3_r by the formula -Lv+(j3+1)*dv

  if(eta[k1] != 0.){
    tmp1_re = (sin(eta[k1]*v1_r) - sin(eta[k1]*v1_l))/eta[k1];
    tmp1_im = (cos(eta[k1]*v1_l) - cos(eta[k1]*v1_r))/eta[k1];
  }
  else 
  {
    tmp1_re = dv;
    tmp1_im = 0.;
  }
  
  if(eta[k2] != 0.){
    tmp2_re = (sin(eta[k2]*v2_r) - sin(eta[k2]*v2_l))/eta[k2];
    tmp2_im = (cos(eta[k2]*v2_l) - cos(eta[k2]*v2_r))/eta[k2];
  }
  else 
  {
    tmp2_re = dv;
    tmp2_im = 0.;
  }
  
  if(eta[k3] != 0.){
    tmp3_re = (sin(eta[k3]*v3_r) - sin(eta[k3]*v3_l))/eta[k3];
    tmp3_im = (cos(eta[k3]*v3_l) - cos(eta[k3]*v3_r))/eta[k3];
  }
  else 
  {
    tmp3_re = dv;
    tmp3_im = 0.;
  }
  
  tem1_re=tmp1_re; tem2_re=tmp2_re;  tem3_re=tmp3_re;  
  tem1_im=tmp1_im; tem2_im=tmp2_im;  tem3_im=tmp3_im;
          
  tmp_re = tmp3_re*(tmp1_re*tmp2_re - tmp1_im*tmp2_im) - tmp3_im*(tmp1_re*tmp2_im + tmp2_re*tmp1_im);
  tmp_im = tmp3_re*(tmp1_re*tmp2_im + tmp2_re*tmp1_im) + tmp3_im*(tmp1_re*tmp2_re - tmp1_im*tmp2_im);
  result[0]=tmp_re; result[1]=tmp_im;
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if(eta[k1] != 0.){
    tmp1_re = ((v1_r*sin(eta[k1]*v1_r) - v1_l*sin(eta[k1]*v1_l))/eta[k1] + (cos(eta[k1]*v1_r) - cos(eta[k1]*v1_l))/eta[k1]/eta[k1] -v1*tem1_re)/dv;
    tmp1_im = ((sin(eta[k1]*v1_r) - sin(eta[k1]*v1_l))/eta[k1]/eta[k1] + (v1_l*cos(eta[k1]*v1_l) - v1_r*cos(eta[k1]*v1_r))/eta[k1] -v1*tem1_im)/dv;
  }
  else
  {
    tmp1_re = 0.; 
    tmp1_im = 0.;
  }         
  
  tmp_re = tem3_re*(tmp1_re*tem2_re - tmp1_im*tem2_im) - tem3_im*(tmp1_re*tem2_im + tem2_re*tmp1_im);
  tmp_im = tem3_re*(tmp1_re*tem2_im + tem2_re*tmp1_im) + tem3_im*(tmp1_re*tem2_re - tmp1_im*tem2_im);
  result[1*2]=tmp_re; result[1*2+1]=tmp_im;
  /////////////////////////////////////////////////////////////////////////////////////////////////////

	  
  if(eta[k2] != 0.){
    tmp2_re = ((v2_r*sin(eta[k2]*v2_r) - v2_l*sin(eta[k2]*v2_l))/eta[k2] + (cos(eta[k2]*v2_r) - cos(eta[k2]*v2_l))/eta[k2]/eta[k2] -v2*tem2_re)/dv;
    tmp2_im = ((sin(eta[k2]*v2_r) - sin(eta[k2]*v2_l))/eta[k2]/eta[k2] + (v2_l*cos(eta[k2]*v2_l) - v2_r*cos(eta[k2]*v2_r))/eta[k2] -v2*tem2_im)/dv;
  }
  else
  {
    tmp2_re = 0; 
    tmp2_im = 0.;
  }

  tmp_re = tem3_re*(tem1_re*tmp2_re - tem1_im*tmp2_im) - tmp3_im*(tem1_re*tmp2_im + tmp2_re*tem1_im);
  tmp_im = tem3_re*(tem1_re*tmp2_im + tmp2_re*tem1_im) + tmp3_im*(tem1_re*tmp2_re - tem1_im*tmp2_im);
  result[2*2]=tmp_re; result[2*2+1]=tmp_im;
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if(eta[k3] != 0.){
    tmp3_re = ((v3_r*sin(eta[k3]*v3_r) - v3_l*sin(eta[k3]*v3_l))/eta[k3] + (cos(eta[k3]*v3_r) - cos(eta[k3]*v3_l))/eta[k3]/eta[k3] -v3*tem3_re)/dv;;
    tmp3_im = ((sin(eta[k3]*v3_r) - sin(eta[k3]*v3_l))/eta[k3]/eta[k3] + (v3_l*cos(eta[k3]*v3_l) - v3_r*cos(eta[k3]*v3_r))/eta[k3] -v3*tem3_im)/dv;;
  }
  else
  {
    tmp3_re = 0.;
    tmp3_im = 0.;
  }
	  
  tmp_re = tmp3_re*(tem1_re*tem2_re - tem1_im*tem2_im) - tmp3_im*(tem1_re*tem2_im + tem2_re*tem1_im);
  tmp_im = tmp3_re*(tem1_re*tem2_im + tem2_re*tem1_im) + tmp3_im*(tem1_re*tem2_re - tem1_im*tem2_im);
  result[3*2]=tmp_re; result[3*2+1]=tmp_im;
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  //BUG: wrote one cos as sin in the following tmp1_im, tmp2_im, tmp3_im !!
  if(eta[k1] != 0.){
    a_re = (v1_r*sin(eta[k1]*v1_r) - v1_l*sin(eta[k1]*v1_l))/eta[k1] + (cos(eta[k1]*v1_r) - cos(eta[k1]*v1_l))/eta[k1]/eta[k1]; // real part of \int exp(i\xi1.v1) v1 dv1
    a_im = (sin(eta[k1]*v1_r) - sin(eta[k1]*v1_l))/eta[k1]/eta[k1] + (v1_l*cos(eta[k1]*v1_l) - v1_r*cos(eta[k1]*v1_r))/eta[k1];
    
    tmp1_re = ((v1_r*v1_r*sin(eta[k1]*v1_r) - v1_l*v1_l*sin(eta[k1]*v1_l) - 2*a_im)/eta[k1] - 2*v1*a_re + v1*v1*tem1_re)/dv/dv;
    tmp1_im = ((v1_l*v1_l*cos(eta[k1]*v1_l) - v1_r*v1_r*cos(eta[k1]*v1_r) + 2*a_re)/eta[k1] - 2*v1*a_im + v1*v1*tem1_im)/dv/dv;
  }
  else
  {
    tmp1_re = dv/12.;
    tmp1_im = 0.;
  }	  
  
  tp1_re = tem3_re*(tmp1_re*tem2_re - tmp1_im*tem2_im) - tem3_im*(tmp1_re*tem2_im + tem2_re*tmp1_im);
  tp1_im = tem3_re*(tmp1_re*tem2_im + tem2_re*tmp1_im) + tem3_im*(tmp1_re*tem2_re - tmp1_im*tem2_im);
  ////////////////////////////////////////////////

  if(eta[k2] != 0.){
    a_re = (v2_r*sin(eta[k2]*v2_r) - v2_l*sin(eta[k2]*v2_l))/eta[k2] + (cos(eta[k2]*v2_r) - cos(eta[k2]*v2_l))/eta[k2]/eta[k2];
    a_im = (sin(eta[k2]*v2_r) - sin(eta[k2]*v2_l))/eta[k2]/eta[k2] + (v2_l*cos(eta[k2]*v2_l) - v2_r*cos(eta[k2]*v2_r))/eta[k2];
    
    tmp2_re = ((v2_r*v2_r*sin(eta[k2]*v2_r) - v2_l*v2_l*sin(eta[k2]*v2_l) - 2*a_im)/eta[k2] - 2*v2*a_re + v2*v2*tem2_re)/dv/dv;
    tmp2_im = ((v2_l*v2_l*cos(eta[k2]*v2_l) - v2_r*v2_r*cos(eta[k2]*v2_r) + 2*a_re)/eta[k2] - 2*v2*a_im + v2*v2*tem2_im)/dv/dv;
  }
  else
  {
    tmp2_re = dv/12.;
    tmp2_im = 0.;
  }         
  
  tp2_re = tem3_re*(tem1_re*tmp2_re - tem1_im*tmp2_im) - tem3_im*(tem1_re*tmp2_im + tmp2_re*tem1_im);
  tp2_im = tem3_re*(tem1_re*tmp2_im + tmp2_re*tem1_im) + tem3_im*(tem1_re*tmp2_re - tem1_im*tmp2_im);
  ////////////////////////////////////////////////  
	  
  if(eta[k3] != 0.){
    a_re = (v3_r*sin(eta[k3]*v3_r) - v3_l*sin(eta[k3]*v3_l))/eta[k3] + (cos(eta[k3]*v3_r) - cos(eta[k3]*v3_l))/eta[k3]/eta[k3];
    a_im = (sin(eta[k3]*v3_r) - sin(eta[k3]*v3_l))/eta[k3]/eta[k3] + (v3_l*cos(eta[k3]*v3_l) - v3_r*cos(eta[k3]*v3_r))/eta[k3];
    
    tmp3_re = ((v3_r*v3_r*sin(eta[k3]*v3_r) - v3_l*v3_l*sin(eta[k3]*v3_l) - 2*a_im)/eta[k3] - 2*v3*a_re + v3*v3*tem3_re)/dv/dv;
    tmp3_im = ((v3_l*v3_l*cos(eta[k3]*v3_l) - v3_r*v3_r*cos(eta[k3]*v3_r) + 2*a_re)/eta[k3] - 2*v3*a_im + v3*v3*tem3_im)/dv/dv;
  }
  else
  {
    tmp3_re = dv/12.;
    tmp3_im = 0.;
  }    
  
  tp3_re = tmp3_re*(tem1_re*tem2_re - tem1_im*tem2_im) - tmp3_im*(tem1_re*tem2_im + tem2_re*tem1_im);
  tp3_im = tmp3_re*(tem1_re*tem2_im + tem2_re*tem1_im) + tmp3_im*(tem1_re*tem2_re - tem1_im*tem2_im);
  ////////////////////////////////////////////////	  
  tmp_re = tp1_re + tp2_re + tp3_re;
  tmp_im = tp1_im + tp2_im + tp3_im;
  result[4*2]=tmp_re; result[4*2+1]=tmp_im;

}

void ProjectedNodeValue(fftw_complex *qHat, double *Q_incremental) // incremental of node values, projected from qHat to DG mesh
{
    int m1,m2,m3,i,j,k,j1,j2,j3, kk,k_v,k_ft, k_eta;
    double tp0, tp2, tp3, tp4, tp5, Q_re, Q_im, u0, u2, u3, u4, u5;
    #pragma omp parallel for private(k_eta, k_ft, k_v, kk, j1,j2,j3,i,j,k,m1,m2,m3,tp0, tp2, tp3, tp4, tp5, Q_re, Q_im, u2, u3, u4, u5) shared(qHat,Q_incremental)
    for(k_ft=0;k_ft<size_ft;k_ft++){
	m3 = k_ft % N; m2 = ((k_ft-m3)/N) % N; m1 = (k_ft - m3 - N*m2)/(N*N);
	j1 = m1*h_v/dv; if(j1==Nv)j1=Nv-1;
	j2 = m2*h_v/dv; if(j2==Nv)j2=Nv-1;
	j3 = m3*h_v/dv; if(j3==Nv)j3=Nv-1; // determine in which element (j1,j2,j3) should this Fourier node (i,j,k) falls
	
	k_v = j1*Nv*Nv + j2*Nv + j3;
	tp0=0.; tp2=0.; tp3=0.; tp4=0.; tp5=0.;
	for(i=0;i<N;i++){
	  for(j=0;j<N;j++){
	    for(k=0;k<N;k++){
	      k_eta = k + N*(j + N*i);
		  //kk = k_v*size_ft+k_eta;
		  IntModes(i,j,k,j1,j2,j3,IntM); //the global IntM must be declared as threadprivate
		  
		  Q_re = nu*qHat[k_eta][0];
	      Q_im = nu*qHat[k_eta][1];
		  //tem = scale3*wtN[i]*wtN[j]*wtN[k]*h_eta*h_eta*h_eta;
		 tp0 += IntM[0]*Q_re - IntM[1]*Q_im;
		  tp2 += IntM[1*2]*Q_re - IntM[1*2+1]*Q_im;
		  tp3 += IntM[2*2]*Q_re - IntM[2*2+1]*Q_im;
		  tp4 += IntM[3*2]*Q_re - IntM[3*2+1]*Q_im;
		  tp5 += IntM[4*2]*Q_re - IntM[4*2+1]*Q_im;
		  
	    }
	  }
	}	      

	u0 = (19*tp0/4. - 15*tp5)/scalev/scaleL/scale3;
	u5 = (60*tp5 - 15*tp0)/scalev/scaleL/scale3;
	u2 = tp2*12./scalev/scaleL/scale3; u3 = tp3*12./scalev/scaleL/scale3; u4 = tp4*12./scalev/scaleL/scale3; 
	
	Q_incremental[k_ft] = u0 + u2*(v[m1]-Gridv((double)j1))/dv + u3*(v[m2]-Gridv((double)j2))/dv + u4*(v[m3]-Gridv((double)j3))/dv + u5*( ((v[m1]-Gridv((double)j1))/dv)*((v[m1]-Gridv((double)j1))/dv) + ((v[m2]-Gridv((double)j2))/dv)*((v[m2]-Gridv((double)j2))/dv) + ((v[m3]-Gridv((double)j3))/dv)*((v[m3]-Gridv((double)j3))/dv) ); 
    }
}	
#ifdef UseMPI
void ComputeQ_FandL(double *f, fftw_complex *qHat, double **conv_weights, fftw_complex *qHat_linear, double **conv_weights_linear)
{
	ComputeQ_MPI_FandL(f, qHat, conv_weights, qHat_linear, conv_weights_linear)
}
void RK4_FandL(double *f, int l, fftw_complex *qHat, double **conv_weights, fftw_complex *qHat_linear, double **conv_weights_linear, double *U, double *dU) //4-th RK. yn=yn+(3*k1+k2+k3+k4)/6
{
	RK4_MPI_FandL(f, l, qHat, conv_weights, qHat_linear, conv_weights_linear, U, dU) //4-th RK. yn=yn+(3*k1+k2+k3+k4)/6
}
void ComputeQ(double *f, fftw_complex *qHat, double **conv_weights)
{
	ComputeQ_MPI(f, qHat, conv_weights)
}
void RK4(double *f, int l, fftw_complex *qHat, double **conv_weights, double *U, double *dU) //4-th RK. yn=yn+(3*k1+k2+k3+k4)/6
{
	RK4_MPI(f, l, qHat, conv_weights, U, dU) //4-th RK. yn=yn+(3*k1+k2+k3+k4)/6
}
#else
void ComputeQ_FandL(double *f, fftw_complex *qHat, double **conv_weights, fftw_complex *qHat_linear, double **conv_weights_linear)
{
	//ComputeQ_NoMPI_FandL(f, qHat, conv_weights, qHat_linear, conv_weights_linear)
}
void RK4_FandL(double *f, int l, fftw_complex *qHat, double **conv_weights, fftw_complex *qHat_linear, double **conv_weights_linear, double *U, double *dU) //4-th RK. yn=yn+(3*k1+k2+k3+k4)/6
{
	//RK4_NoMPI_FandL(f, l, qHat, conv_weights, qHat_linear, conv_weights_linear, U) //4-th RK. yn=yn+(3*k1+k2+k3+k4)/6
}
void ComputeQ(double *f, fftw_complex *qHat, double **conv_weights)
{
	ComputeQ_NoMPI(f, qHat, conv_weights)
}
void RK4(double *f, int l, fftw_complex *qHat, double **conv_weights, double *U, double *dU) //4-th RK. yn=yn+(3*k1+k2+k3+k4)/6
{
	RK4_NoMPI(f, l, qHat, conv_weights, U) //4-th RK. yn=yn+(3*k1+k2+k3+k4)/6
}

//#ifdef UseMPI
//#ifdef FullandLinear
void ComputeQ_MPI_FandL(double *f, fftw_complex *qHat, double **conv_weights, fftw_complex *qHat_linear, double **conv_weights_linear)
{
  int i, j, k, l, m, n, x, y, z;
  int start_i, start_j, start_k, end_i, end_j, end_k;
  //fftw_complex *fftIn, *fftOut;
  double tempD, tempD1, tmp0, tmp1, tmp01, tmp11;
  double prefactor = h_eta*h_eta*h_eta; // we dont have scale3 here.

  //fftIn = (fftw_complex *)fftw_malloc(N*N*N*sizeof(fftw_complex));
  //fftOut = (fftw_complex *)fftw_malloc(N*N*N*sizeof(fftw_complex));
  
  for(i=0;i<size_ft;i++){
    fftIn[i][0] = f[i];
    fftIn[i][1] = 0.;
  }

  fft3D(fftIn, fftOut);
  
  //printf("fft done\n");
  #pragma omp parallel for private(j,k,l,m,n,x,y,z,start_i,start_j,start_k,end_i,end_j,end_k,tempD, tempD1,tmp0, tmp1) shared(qHat, qHat_linear, fftOut, conv_weights, conv_weights_linear)
  for(i=0;i<N;i++) {
    for(j=0;j<N;j++){
      for(k=0;k<N;k++)
	{
	  //figure out the windows for the convolutions (i.e. where xi(l) and eta(i)-xi(l) are in the domain)
	  if( i < N/2 ) {
	    start_i = 0;
	    end_i = i + N/2 + 1; 
	  }
	  else {
	    start_i = i - N/2 + 1;
	    end_i = N;
	  }
	  
	  if( j < N/2 ) {
	    start_j = 0;
	    end_j = j + N/2 + 1; 
	  }
	  else {
	    start_j = j - N/2 + 1;
	    end_j = N;
	  }
	  
	  if( k < N/2 ) {
	    start_k = 0;
	    end_k = k + N/2 + 1; 
	  }
	  else {
	    start_k = k - N/2 + 1;
	    end_k = N;
	  }
	  tmp0=0.; tmp1=0.; tmp01=0.; tmp11=0.;
	  for(l=start_i;l<end_i;l++) {
	    for(m=start_j;m<end_j;m++) {
	      for(n=start_k;n<end_k;n++)
		{
		
		  x = i + N/2 - l;
		  y = j + N/2 - m;
		  z = k + N/2 - n;
		  //printf("%d, %d, %d; %d, %d, %d\n",i,j,k,l,m,n);
		 // printf("tempD=%g\n",conv_weights[k + N*(j+ N*i)][n + N*(m + N*l)]);
		  //get convolution weight		   
		  tempD = conv_weights[k + N*(j+ N*i)][n + N*(m + N*l)];
		  tempD1 = conv_weights_linear[k + N*(j+ N*i)][n + N*(m + N*l)];		  
		//  printf(" tempD=%g, prefactor=%g, fftOut=[%g,%g] at %d, %d, %d; %d, %d, %d; %d, %d, %d\n",tempD, prefactor, fftOut[z + N*(y + N*x)][0], fftOut[z + N*(y + N*x)][1],i,j,k,l,m,n, x,y,z);
		  tmp0 += prefactor*wtN[l]*wtN[m]*wtN[n]*(tempD*(fftOut[n + N*(m + N*l)][0]*fftOut[z + N*(y + N*x)][0] - fftOut[n + N*(m + N*l)][1]*fftOut[z + N*(y + N*x)][1]) + scale3*tempD1*fftOut[z + N*(y + N*x)][0]);
		  tmp1 += prefactor*wtN[l]*wtN[m]*wtN[n]*(tempD*(fftOut[n + N*(m + N*l)][0]*fftOut[z + N*(y + N*x)][1] + fftOut[n + N*(m + N*l)][1]*fftOut[z + N*(y + N*x)][0]) + scale3*tempD1*fftOut[z + N*(y + N*x)][1]);
		  
		  tmp01 += prefactor*wtN[l]*wtN[m]*wtN[n]*scale3*tempD1*fftOut[z + N*(y + N*x)][0];
		  tmp11 += prefactor*wtN[l]*wtN[m]*wtN[n]*scale3*tempD1*fftOut[z + N*(y + N*x)][1];
		 
		}
	    }
	  }
	// printf("%d, %d, %d done\n", i,j,k);
	 qHat[k + N*(j + N*i)][0] = tmp0;  
	 qHat[k + N*(j + N*i)][1] = tmp1;
	 qHat_linear[k + N*(j + N*i)][0] = tmp01;  
	 qHat_linear[k + N*(j + N*i)][1] = tmp11;
	}
    }
  }

}

void RK4_MPI_FandL(double *f, int l, fftw_complex *qHat, double **conv_weights, fftw_complex *qHat_linear, double **conv_weights_linear, double *U, double *dU) //4-th RK. yn=yn+(3*k1+k2+k3+k4)/6
{
int i,j,k, j1, j2, j3, k_v, k_eta, kk, l_local;  
  double Q_re, Q_im, tp0, tp2, tp3,tp4,tp5, tmp0=0., tmp2=0., tmp3=0., tmp4=0.,tmp5=0., tem;

  l_local = l%chunk_Nx;
  
  #pragma omp parallel for private(i) shared(qHat, qHat_linear)
  for(i=0;i<size_ft;i++){
    qHat[i][0] += qHat_linear[i][0];
	qHat[i][1] += qHat_linear[i][1];
  }
  FS(qHat, fftOut); 

  #pragma omp parallel for private(i) shared(Q,fftOut,f1,f)
  for(i=0;i<size_ft;i++){    
    Q[i] = fftOut[i][0];
    f1[i] = f[i] + dt*Q[i]*nu; //BUG: this evolution (only on node values) is not consistent with our conservation routine, which preserves the exact moments of the {1,v,|v|^{2}} approximations
  }

  ComputeQ_FandL(f1, Q1_fft, conv_weights, Q1_fft_linear, conv_weights_linear);
  conserveAllMoments(Q1_fft, Q1_fft_linear);   	//conserves k2	
  
  #pragma omp parallel for private(i) shared(Q1_fft, Q1_fft_linear)
  for(i=0;i<size_ft;i++){
    Q1_fft[i][0] += Q1_fft_linear[i][0];
	Q1_fft[i][1] += Q1_fft_linear[i][1];
  }
  
  FS(Q1_fft, fftOut);

  #pragma omp parallel for private(i) shared(Q,Q1,fftOut,f1,f)
  for(i=0;i<size_ft;i++){ 	
    Q1[i] = fftOut[i][0];	   
    f1[i] = f[i] +  0.5*dt*Q[i]*nu + 0.5*dt*Q1[i]*nu;
  }
 
  ComputeQ_FandL(f1, Q2_fft, conv_weights, Q2_fft_linear, conv_weights_linear);
  conserveAllMoments(Q2_fft, Q2_fft_linear);   //conserves k3

  #pragma omp parallel for private(i) shared(Q2_fft, Q2_fft_linear)
  for(i=0;i<size_ft;i++){
    Q2_fft[i][0] += Q2_fft_linear[i][0];
	Q2_fft[i][1] += Q2_fft_linear[i][1];
  }
  FS(Q2_fft, fftOut);
  //ifft3D(Q2_fft, fftOut);
  #pragma omp parallel for private(i) shared(Q,Q1,fftOut,f1,f)
  for(i=0;i<size_ft;i++){	
    Q1[i] = fftOut[i][0];	 
    f1[i] = f[i] + 0.5*dt*Q[i]*nu + 0.5*dt*Q1[i]*nu;
  }

  ComputeQ_FandL(f1, Q3_fft, conv_weights, Q3_fft_linear, conv_weights_linear);
  conserveAllMoments(Q3_fft, Q3_fft_linear);                //conserves k4

  #pragma omp parallel for private(i) shared(Q3_fft, Q3_fft_linear)
  for(i=0;i<size_ft;i++){
    Q3_fft[i][0] += Q3_fft_linear[i][0];
	Q3_fft[i][1] += Q3_fft_linear[i][1];
  }
  #pragma omp parallel for schedule(dynamic) private(j1,j2,j3,i,j,k,k_v,k_eta,Q_re, Q_im, tp0, tp2,tp3,tp4, tp5) shared(l,l_local,qHat,U, dU) //reduction(+: tmp0, tmp2, tmp3, tmp4, tmp5)
  for(int kt=0;kt<size_v;kt++){
    j3 = kt % Nv; j2 = ((kt-j3)/Nv) % Nv; j1 = (kt - j3 - Nv*j2)/(Nv*Nv);
    tp0=0.; tp2=0.; tp3=0.; tp4=0.; tp5=0.;
    for(i=0;i<N;i++){
      for(j=0;j<N;j++){
		for(k=0;k<N;k++){
		  k_eta = k + N*(j + N*i);

		  IntModes(i,j,k,j1,j2,j3,IntM); //the global IntM must be declared as threadprivate; BUG: forgot to uncomment this and thus IntM=0 !!
		  
		  Q_re = nu*(0.5*qHat[k_eta][0] + (Q1_fft[k_eta][0]+Q2_fft[k_eta][0]+Q3_fft[k_eta][0])/6.);
		  Q_im = nu*(0.5*qHat[k_eta][1] + (Q1_fft[k_eta][1]+Q2_fft[k_eta][1]+Q3_fft[k_eta][1])/6.);
		  
		  //tem = scale3*wtN[i]*wtN[j]*wtN[k]*h_eta*h_eta*h_eta;
		  tp0 += IntM[0]*Q_re - IntM[1]*Q_im;
		  tp2 += IntM[1*2]*Q_re - IntM[1*2+1]*Q_im;
		  tp3 += IntM[2*2]*Q_re - IntM[2*2+1]*Q_im;
		  tp4 += IntM[3*2]*Q_re - IntM[3*2+1]*Q_im;
		  tp5 += IntM[4*2]*Q_re - IntM[4*2+1]*Q_im;
		}
      }
    }	     
    //tmp0 += tp0; tmp2 += dv*tp2 + Gridv((double)j1)*tp0; tmp3 += dv*tp3 +Gridv((double)j2)*tp0 ;tmp4 += dv*tp4 + Gridv((double)j3)*tp0;  //tmp5 += (Gridv((double)j1)*Gridv((double)j1) + Gridv((double)j2)*Gridv((double)j2) + Gridv((double)j3)*Gridv((double)j3))*tp0 +dv*dv*tp5+2*dv*(Gridv((double)j1)*tp2 + Gridv((double)j2)*tp3 +Gridv((double)j3)*tp4);     
	//tmp5 += dv*dv*tp5 -  (Gridv((double)j1)*Gridv((double)j1) + Gridv((double)j2)*Gridv((double)j2) + Gridv((double)j3)*Gridv((double)j3))*tp0 + 2*(Gridv((double)j1)*(dv*tp2 + Gridv((double)j1)*tp0) + Gridv((double)j2)*(dv*tp3 + Gridv((double)j2)*tp0) + Gridv((double)j3)*(dv*tp4 + Gridv((double)j3)*tp0));
	  
    k_v = l*size_v + kt;      
    tp0 = U[k_v*6+0] + U[k_v*6+5]/4. + dt*tp0/scalev/scaleL/scale3;
    tp2 = U[k_v*6+2] + dt*tp2*12./scalev/scaleL/scale3;
    tp3 = U[k_v*6+3] + dt*tp3*12./scalev/scaleL/scale3;
    tp4 = U[k_v*6+4] + dt*tp4*12./scalev/scaleL/scale3;
    tp5 = U[k_v*6+0]/4. + U[k_v*6+5]*19./240. + dt*tp5/scalev/scaleL/scale3;

    dU[(l_local*size_v + kt)*5] = 19*tp0/4. - 15*tp5;
    dU[(l_local*size_v + kt)*5+4] = 60*tp5 - 15*tp0;
    dU[(l_local*size_v + kt)*5+1] = tp2;
    dU[(l_local*size_v + kt)*5+2] = tp3;
    dU[(l_local*size_v + kt)*5+3] = tp4;   
  }
}
//#else
void ComputeQ_MPI(double *f, fftw_complex *qHat, double **conv_weights)
{
	int i, j, k, l, m, n, x, y, z;												// declare (i,j,k) (the indices for a given value of given ki = ki_(i,j,k)), (l,m,n) (counters for the quadrature to calculate the integral w.r.t. eta in the evaluation of qHat and so also represent the indices of a given eta = eta_(l,m,n)) & (x,y,z) (the indices for the value of a subtraction in the calculation, namely eta_(x,y,z) = ki_(i,j,k) - eta_(l,m,n))
	int start_i, start_j, start_k, end_i, end_j, end_k;							// declare start_i, start_j & start_k (the indices for the values of the lower bounds of integration in computation of the convolution, corresponding to the lowest point where both functions are non-zero, in each velocity direction) and end_i, end_j & end_k (the indices for the values of the upper bounds of integration in computation of the convolution, corresponding to the highest point where both functions are non-zero, in each velocity direction)
	double tempD, tmp0, tmp1;													// declare tempD (the value of the convolution weight at a given ki & eta), tmp0 (which will become the real part of qHat) & tmp1 (which will become the imaginary part of qHat)
	double prefactor = h_eta*h_eta*h_eta; 										// declare prefactor (the value of h_eta^3, as no scale3 in Fourier space) and set its value
  
	for(i=0;i<size_ft;i++)														// initialise the input of the FFT
	{
		fftIn[i][0] = f[i];														// set the real part to the sampling of the solution stored in f
		fftIn[i][1] = 0.;														// set the imaginary part to zero
	}

	fft3D(fftIn, fftOut);														// perform the FFT of fftIn and store the result in fftOut
  
	//printf("fft done\n");
	#pragma omp parallel for schedule(dynamic) private(i,j,k,l,m,n,x,y,z,start_i,start_j,start_k,end_i,end_j,end_k,tempD) shared(qHat, fftOut, conv_weights) reduction(+:tmp0, tmp1)
	for(i=0;i<N;i++) 															// loop through all points in the ki_1
	{
		for(j=0;j<N;j++)														// loop through all points in the ki_2
		{
			for(k=0;k<N;k++)													// loop through all points in the ki_3
			{
				//figure out the windows for the convolutions (i.e. where eta(l,m,n) and ki(i,j,k)-eta(l,m,n) are in the domain)
				if( i < N/2 ) 													// if ki_1(i) < 0 then the values of l for which f(ki_1(i)-eta_1(l))*f(eta_1(l)) give a non-zero product range need -Lv < eta_1 <= ki_1 + Lv/2, due to the support of f being -Lv to Lv
				{
					start_i = 0;												// set start_i to 0 to represent -Lv as the lower bound of integration
					end_i = i + N/2 + 1;										// set end_i to i + N/2 + 1 to represent eta_1(i+N/2+1/2) as the upper bound of integration (with +1 there so that all cells less than index end_i are integrated over)
				}
				else 															// if ki_1(i) >= 0 then the values of l for which f(ki_1(i)-eta_1(l))*f(eta_1(l)) give a non-zero product range need ki_1 - Lv/2 < eta_1 <= Lv, due to the support of f being -Lv to Lv
				{
					start_i = i - N/2 + 1;										// set start_i to i - N/2 + 1 to represent eta_1(i) - Lv as the lower bound of integration (with +1 since ki_1[i-N/2] - eta_1[i] actually overlaps with eta_1[-1/2] = -Lv, where f(-Lv)=0, so start at the next index)
					end_i = N;													// set end_i to N to represent Lv as the upper bound of integration (i.e. integrate over all cells with index less than N)
				}
	  
				if( j < N/2 )													// if ki_2(j) < 0 then the values of m for which f(ki_2(j)-eta_2(m))*f(eta_2(m)) give a non-zero product range need -Lv < eta_2 <= ki_2 + Lv/2, due to the support of f being -Lv to Lv
				{
					start_j = 0;												// set start_j to 0 to represent -Lv as the lower bound of integration
					end_j = j + N/2 + 1;										// set end_j to j + N/2 + 1 to represent eta_2(j+N/2+1/2) as the upper bound of integration (with +1 there so that all cells less than index end_i are integrated over)
				}
				else															// if ki_2(j) >= 0 then the values of m for which f(ki_2(j)-eta_2(m))*f(eta_2(m)) give a non-zero product range need ki_2 - Lv/2 < eta_2 <= Lv, due to the support of f being -Lv to Lv
				{
					start_j = j - N/2 + 1;										// set start_j to j - N/2 + 1 to represent eta_2(j) - Lv as the lower bound of integration (with +1 since ki_2[j-N/2] - eta_2[j] actually overlaps with eta_2[-1/2] = -Lv, where f(-Lv)=0, so start at the next index)
					end_j = N;													// set end_j to N to represent Lv as the upper bound of integration (i.e. integrate over all cells with index less than N)
				}
	  
				if( k < N/2 )													// if ki_3(k) < 0 then the values of n for which f(ki_3(k)-eta_3(n))*f(eta_3(n)) give a non-zero product range need -Lv < eta_3 <= ki_3 + Lv/2, due to the support of f being -Lv to Lv
				{
					start_k = 0;												// set start_k to 0 to represent -Lv as the lower bound of integration
					end_k = k + N/2 + 1;										// set end_k to k + N/2 + 1 to represent eta_3(k+N/2+1/2) as the upper bound of integration (with +1 there so that all cells less than index end_i are integrated over)
				}
				else															// if ki_3(k) >= 0 then the values of n for which f(ki_3(k)-eta_3(n))*f(eta_3(n)) give a non-zero product range need ki_3 - Lv/2 < eta_3 <= Lv, due to the support of f being -Lv to Lv
				{
					start_k = k - N/2 + 1;										// set start_k to k - N/2 + 1 to represent eta_3(k) - Lv as the lower bound of integration (with +1 since ki_3[k-N/2] - eta_3[k] actually overlaps with eta_3[-1/2] = -Lv, where f(-Lv)=0, so start at the next index)
					end_k = N;													// set end_k to N to represent Lv as the upper bound of integration (i.e. integrate over all cells with index less than N)
				}
				tmp0=0.; tmp1=0.;												// initialise tmp0 & tmp1 at zero to begin the quadrature
				for(l=start_i;l<end_i;l++)										// loop through all the quadrature indices in the eta_1 direction that give non-zero contribution
				{
					for(m=start_j;m<end_j;m++)									// loop through all the quadrature indices in the eta_2 direction that give non-zero contribution
					{
						for(n=start_k;n<end_k;n++)								// loop through all the quadrature indices in the eta_3 direction that give non-zero contribution
						{
		
							x = i + N/2 - l;									// set the index x to i + N/2 - l to represent the subtraction eta[x] = ki[i] - eta[l]
							y = j + N/2 - m;									// set the index y to j + N/2 - m to represent the subtraction eta[y] = ki[j] - eta[m]
							z = k + N/2 - n;									// set the index z to k + N/2 - n to represent the subtraction eta[z] = ki[k] - eta[n]
							//printf("%d, %d, %d; %d, %d, %d\n",i,j,k,l,m,n);
							// printf("tempD=%g\n",conv_weights[k + N*(j+ N*i)][n + N*(m + N*l)]);
							tempD = conv_weights[k + N*(j+ N*i)][n + N*(m + N*l)];		// set tempD to the value of the convolution weight corresponding to current value of ki(i,j,k) & eta(l,m,n)
							tmp0 += prefactor*wtN[l]*wtN[m]*wtN[n]*tempD*(fftOut[n + N*(m + N*l)][0]*fftOut[z + N*(y + N*x)][0] - fftOut[n + N*(m + N*l)][1]*fftOut[z + N*(y + N*x)][1]);		// increment the value of the real part of qHat(ki(i,j,k)) by fHat(eta(l,m,n))*f(ki(i,j,k)-eta(l,m,n))*conv_weight(ki(i,j,k),eta(l,m,n)) for the current values of l, m & n in the quadrature sum
		  
							tmp1 += prefactor*wtN[l]*wtN[m]*wtN[n]*tempD*(fftOut[n + N*(m + N*l)][0]*fftOut[z + N*(y + N*x)][1] + fftOut[n + N*(m + N*l)][1]*fftOut[z + N*(y + N*x)][0]);		// increment the value of the imaginary part of qHat(ki(i,j,k)) by fHat(eta(l,m,n))*f(ki(i,j,k)-eta(l,m,n))*conv_weight(ki(i,j,k),eta(l,m,n)) for the current values of l, m & n in the quadrature sum
							// printf(" tempD=%g, prefactor=%g, fftOut=[%g,%g] at %d, %d, %d; %d, %d, %d; %d, %d, %d\n",tempD, prefactor, fftOut[z + N*(y + N*x)][0], fftOut[z + N*(y + N*x)][1],i,j,k,l,m,n, x,y,z);
						}
					}
				}
				// printf("%d, %d, %d done\n", i,j,k);
				qHat[k + N*(j + N*i)][0] = tmp0;								// set the real part of qHat(ki(i,j,k)) to the value tmp0 calculated in the quadrature
				qHat[k + N*(j + N*i)][1] = tmp1;								// set the imaginary part of qHat(ki(i,j,k)) to the value tmp1 calculated in the quadrature
				// printf("%d, %d, %d write-in done\n", i,j,k);
			}
		}
	}
}

void RK4_MPI(double *f, int l, fftw_complex *qHat, double **conv_weights, double *U, double *dU) //4-th RK. yn=yn+(3*k1+k2+k3+k4)/6
{
  int i,j,k, j1, j2, j3, k_v, k_eta, kk,l_local;  
  double Q_re, Q_im, tp0, tp2, tp3,tp4,tp5, tmp0=0., tmp2=0., tmp3=0., tmp4=0.,tmp5=0., tem;

  l_local = l%chunk_Nx;
  
  FS(qHat, fftOut); 																	// set fftOut to the Fourier series representation of qHat (i.e. the IFFT of qHat)
  //ifft3D(qHat, fftOut);
  #pragma omp parallel for private(i) shared(Q,fftOut,f1,f)
  for(i=0;i<size_ft;i++)																// calculate the first step of RK4
  {
    Q[i] = fftOut[i][0];																// this is Q(Fn, Fn) so that Kn^1 = dt*Q(Fn, Fn) = dt*Q[i]
    f1[i] = f[i] + dt*Q[i]*nu; 															// this is Fn + Kn^1(*nu...?) BUG: this evolution (only on node values) is not consistent with our conservation routine, which preserves the exact moments of the {1,v,|v|^{2}} approximations
  }

  ComputeQ(f1, Q1_fft, conv_weights);													// calculate the Fourier tranform of Q(f1,f1) using conv_weights for the weights in the convolution, then store the results of the Fourier transform in Q1_fft
  conserveAllMoments(Q1_fft);   														// perform the explicit conservation calculation on Kn2^ = Q^(f1,f1) = Q1_fft

  FS(Q1_fft, fftOut);																	// set fftOut to the Fourier series representation of Q1_fft (i.e. the IFFT of Q1_fft, so that Kn^2 = fftOut = Q(Fn + dt*Kn^1, Fn + dt*Kn^1) )
  //ifft3D(Q1_fft, fftOut);
  #pragma omp parallel for private(i) shared(Q,Q1,fftOut,f1,f)
  for(i=0;i<size_ft;i++)																// calculate the second step of RK4
  {
    Q1[i] = fftOut[i][0];	   
    f1[i] = f[i] +  0.5*dt*Q[i]*nu + 0.5*dt*Q1[i]*nu;
  }

  ComputeQ(f1, Q2_fft, conv_weights); //collides
  conserveAllMoments(Q2_fft);   //conserves k3

  FS(Q2_fft, fftOut);
  //ifft3D(Q2_fft, fftOut);
  #pragma omp parallel for private(i) shared(Q,Q1,fftOut,f1,f)
  for(i=0;i<size_ft;i++)																// calculate the third step of RK4
  {
    Q1[i] = fftOut[i][0];	 
    f1[i] = f[i] + 0.5*Q[i]*nu + 0.5*Q1[i]*nu;
  }

  ComputeQ(f1, Q3_fft, conv_weights); //collides
  conserveAllMoments(Q3_fft);                //conserves k4

  #pragma omp parallel for schedule(dynamic) private(j1,j2,j3,i,j,k,k_v,k_eta,kk,Q_re, Q_im) shared(l, l_local, qHat,U, dU) reduction(+:tp0, tp2,tp3,tp4, tp5)  // calculate the fourth step of RK4 (still in Fourier space though?!) - reduction(+: tmp0, tmp2, tmp3, tmp4, tmp5)
  for(int kt=0;kt<size_v;kt++){
    j3 = kt % Nv; j2 = ((kt-j3)/Nv) % Nv; j1 = (kt - j3 - Nv*j2)/(Nv*Nv);
    tp0=0.; tp2=0.; tp3=0.; tp4=0.; tp5=0.;
    for(i=0;i<N;i++){
      for(j=0;j<N;j++){
		for(k=0;k<N;k++){
		  k_eta = k + N*(j + N*i);
		  
		  IntModes(i,j,k,j1,j2,j3,IntM); //the global IntM must be declared as threadprivate; BUG: forgot to uncomment this and thus IntM=0 !!
		  
		  Q_re = nu*(0.5*qHat[k_eta][0] + (Q1_fft[k_eta][0]+Q2_fft[k_eta][0]+Q3_fft[k_eta][0])/6.);
		  Q_im = nu*(0.5*qHat[k_eta][1] + (Q1_fft[k_eta][1]+Q2_fft[k_eta][1]+Q3_fft[k_eta][1])/6.);
		  
		  //tem = scale3*wtN[i]*wtN[j]*wtN[k]*h_eta*h_eta*h_eta;
		  tp0 += IntM[0]*Q_re - IntM[1]*Q_im;
		  tp2 += IntM[1*2]*Q_re - IntM[1*2+1]*Q_im;
		  tp3 += IntM[2*2]*Q_re - IntM[2*2+1]*Q_im;
		  tp4 += IntM[3*2]*Q_re - IntM[3*2+1]*Q_im;
		  tp5 += IntM[4*2]*Q_re - IntM[4*2+1]*Q_im;
		}
      }
    }	     
    //tmp0 += tp0; tmp2 += dv*tp2 + Gridv((double)j1)*tp0; tmp3 += dv*tp3 +Gridv((double)j2)*tp0 ;tmp4 += dv*tp4 + Gridv((double)j3)*tp0;  //tmp5 += (Gridv((double)j1)*Gridv((double)j1) + Gridv((double)j2)*Gridv((double)j2) + Gridv((double)j3)*Gridv((double)j3))*tp0 +dv*dv*tp5+2*dv*(Gridv((double)j1)*tp2 + Gridv((double)j2)*tp3 +Gridv((double)j3)*tp4);     
	//tmp5 += dv*dv*tp5 -  (Gridv((double)j1)*Gridv((double)j1) + Gridv((double)j2)*Gridv((double)j2) + Gridv((double)j3)*Gridv((double)j3))*tp0 + 2*(Gridv((double)j1)*(dv*tp2 + Gridv((double)j1)*tp0) + Gridv((double)j2)*(dv*tp3 + Gridv((double)j2)*tp0) + Gridv((double)j3)*(dv*tp4 + Gridv((double)j3)*tp0));
	  
    k_v = l*size_v + kt;      
    tp0 = U[k_v*6+0] + U[k_v*6+5]/4. + dt*tp0/scalev/scaleL/scale3;
    tp2 = U[k_v*6+2] + dt*tp2*12./scalev/scaleL/scale3;
    tp3 = U[k_v*6+3] + dt*tp3*12./scalev/scaleL/scale3;
    tp4 = U[k_v*6+4] + dt*tp4*12./scalev/scaleL/scale3;
    tp5 = U[k_v*6+0]/4. + U[k_v*6+5]*19./240. + dt*tp5/scalev/scaleL/scale3;

    dU[(l_local*size_v + kt)*5] = 19*tp0/4. - 15*tp5;
    dU[(l_local*size_v + kt)*5+4] = 60*tp5 - 15*tp0;
    dU[(l_local*size_v + kt)*5+1] = tp2;
    dU[(l_local*size_v + kt)*5+2] = tp3;
    dU[(l_local*size_v + kt)*5+3] = tp4;   
  }
}
//#endif
 
//#else

//#ifdef FullandLinear
/*void ComputeQ_NoMPI_FandL(double *f, fftw_complex *qHat, double **conv_weights, fftw_complex *qHat_linear, double **conv_weights_linear)
{
  int i, j, k, l, m, n, x, y, z;
  int start_i, start_j, start_k, end_i, end_j, end_k;
  //fftw_complex *fftIn, *fftOut;
  double tempD, tempD1, tmp0, tmp1;
  double prefactor = h_eta*h_eta*h_eta; // we dont have scale3 here.

  //fftIn = (fftw_complex *)fftw_malloc(N*N*N*sizeof(fftw_complex));
  //fftOut = (fftw_complex *)fftw_malloc(N*N*N*sizeof(fftw_complex));
  
  for(i=0;i<size_ft;i++){
    fftIn[i][0] = f[i];
    fftIn[i][1] = 0.;
  }

  fft3D(fftIn, fftOut);
  
  //printf("fft done\n");
  #pragma omp parallel for private(j,k,l,m,n,x,y,z,start_i,start_j,start_k,end_i,end_j,end_k,tempD, tempD1,tmp0, tmp1) shared(qHat, fftOut, conv_weights, conv_weights_linear)
  for(i=0;i<N;i++) {
    for(j=0;j<N;j++){
      for(k=0;k<N;k++)
	{
	  //figure out the windows for the convolutions (i.e. where xi(l) and eta(i)-xi(l) are in the domain)
	  if( i < N/2 ) {
	    start_i = 0;
	    end_i = i + N/2 + 1; 
	  }
	  else {
	    start_i = i - N/2 + 1;
	    end_i = N;
	  }
	  
	  if( j < N/2 ) {
	    start_j = 0;
	    end_j = j + N/2 + 1; 
	  }
	  else {
	    start_j = j - N/2 + 1;
	    end_j = N;
	  }
	  
	  if( k < N/2 ) {
	    start_k = 0;
	    end_k = k + N/2 + 1; 
	  }
	  else {
	    start_k = k - N/2 + 1;
	    end_k = N;
	  }
	  tmp0=0.; tmp1=0.;
	  for(l=start_i;l<end_i;l++) {
	    for(m=start_j;m<end_j;m++) {
	      for(n=start_k;n<end_k;n++)
		{
		
		  x = i + N/2 - l;
		  y = j + N/2 - m;
		  z = k + N/2 - n;
		  //printf("%d, %d, %d; %d, %d, %d\n",i,j,k,l,m,n);
		 // printf("tempD=%g\n",conv_weights[k + N*(j+ N*i)][n + N*(m + N*l)]);
		  //get convolution weight		   
		  tempD = conv_weights[k + N*(j+ N*i)][n + N*(m + N*l)];
		  tempD1 = conv_weights_linear[k + N*(j+ N*i)][n + N*(m + N*l)];		  
		//  printf(" tempD=%g, prefactor=%g, fftOut=[%g,%g] at %d, %d, %d; %d, %d, %d; %d, %d, %d\n",tempD, prefactor, fftOut[z + N*(y + N*x)][0], fftOut[z + N*(y + N*x)][1],i,j,k,l,m,n, x,y,z);
		  tmp0 += prefactor*wtN[l]*wtN[m]*wtN[n]*(tempD*(fftOut[n + N*(m + N*l)][0]*fftOut[z + N*(y + N*x)][0] - fftOut[n + N*(m + N*l)][1]*fftOut[z + N*(y + N*x)][1]) + scale3*tempD1*fftOut[z + N*(y + N*x)][0]);
		  
		  tmp1 += prefactor*wtN[l]*wtN[m]*wtN[n]*(tempD*(fftOut[n + N*(m + N*l)][0]*fftOut[z + N*(y + N*x)][1] + fftOut[n + N*(m + N*l)][1]*fftOut[z + N*(y + N*x)][0]) + scale3*tempD1*fftOut[z + N*(y + N*x)][1]);
		 // printf(" tempD=%g, prefactor=%g, fftOut=[%g,%g] at %d, %d, %d; %d, %d, %d; %d, %d, %d\n",tempD, prefactor, fftOut[z + N*(y + N*x)][0], fftOut[z + N*(y + N*x)][1],i,j,k,l,m,n, x,y,z);
		}
	    }
	  }
	// printf("%d, %d, %d done\n", i,j,k);
	 qHat[k + N*(j + N*i)][0] = tmp0;  
	 qHat[k + N*(j + N*i)][1] = tmp1;
	 // printf("%d, %d, %d write-in done\n", i,j,k);
	}
    }
  }

 // fftw_free(fftIn);
 // fftw_free(fftOut);

}
void RK4_NoMPI_FandL(double *f, int l, fftw_complex *qHat, double **conv_weights, fftw_complex *qHat_linear, double **conv_weights_linear, double *U) //4-th RK. yn=yn+(3*k1+k2+k3+k4)/6
{
  int i,j,k, j1, j2, j3, k_v, k_eta,kk;  
  double Q_re, Q_im, tp0, tp2, tp3,tp4,tp5, tmp0=0., tmp2=0., tmp3=0., tmp4=0.,tmp5=0., tem;

  FS(qHat, fftOut); 
  //ifft3D(qHat, fftOut);
  #pragma omp parallel for private(i) shared(Q,fftOut,f1,f)
  for(i=0;i<size_ft;i++){    
    Q[i] = fftOut[i][0];
    f1[i] = f[i] + dt*Q[i]*nu; //BUG: this evolution (only on node values) is not consistent with our conservation routine, which preserves the exact moments of the {1,v,|v|^{2}} approximations
  }

  ComputeQ_FandL(f1, Q1_fft, conv_weights, qHat_linear, conv_weights_linear); //collides
  conserveAllMoments(Q1_fft);   	//conserves k2	

  FS(Q1_fft, fftOut);
  //ifft3D(Q1_fft, fftOut);
  #pragma omp parallel for private(i) shared(Q,Q1,fftOut,f1,f)
  for(i=0;i<size_ft;i++){ 	
    Q1[i] = fftOut[i][0];	   
    f1[i] = f[i] +  0.5*dt*Q[i]*nu + 0.5*dt*Q1[i]*nu;
  }
  
  ComputeQ_FandL(f1, Q2_fft, conv_weights, qHat_linear, conv_weights_linear); //collides
  conserveAllMoments(Q2_fft);   //conserves k3

  FS(Q2_fft, fftOut);
  //ifft3D(Q2_fft, fftOut);
  #pragma omp parallel for private(i) shared(Q,Q1,fftOut,f1,f)
  for(i=0;i<size_ft;i++){	
    Q1[i] = fftOut[i][0];	 
    f1[i] = f[i] + 0.5*dt*Q[i]*nu + 0.5*dt*Q1[i]*nu;
  }
	 
  ComputeQ_FandL(f1, Q3_fft, conv_weights, qHat_linear, conv_weights_linear); //collides
  conserveAllMoments(Q3_fft);                //conserves k4
  
  #pragma omp parallel for schedule(dynamic) private(j1,j2,j3,i,j,k,k_v,k_eta,kk,Q_re, Q_im, tp0, tp2,tp3,tp4, tp5) shared(l,qHat,U) //reduction(+: tmp0, tmp2, tmp3, tmp4, tmp5)
  for(int kt=0;kt<size_v;kt++){
    j3 = kt % Nv; j2 = ((kt-j3)/Nv) % Nv; j1 = (kt - j3 - Nv*j2)/(Nv*Nv);
    tp0=0.; tp2=0.; tp3=0.; tp4=0.; tp5=0.;
    for(i=0;i<N;i++){
      for(j=0;j<N;j++){
		for(k=0;k<N;k++){
		  k_eta = k + N*(j + N*i);
		  //kk = kt*size_ft+k_eta;
		  //IntModes(i,j,k,j1,j2,j3,IntM); //the global IntM must be declared as threadprivate
		  
		  Q_re = nu*(0.5*qHat[k_eta][0] + (Q1_fft[k_eta][0]+Q2_fft[k_eta][0]+Q3_fft[k_eta][0])/6.);
		  Q_im = nu*(0.5*qHat[k_eta][1] + (Q1_fft[k_eta][1]+Q2_fft[k_eta][1]+Q3_fft[k_eta][1])/6.);
		  
		  //tem = scale3*wtN[i]*wtN[j]*wtN[k]*h_eta*h_eta*h_eta;
		  tp0 += IntM[0]*Q_re - IntM[1]*Q_im;
		  tp2 += IntM[1*2]*Q_re - IntM[1*2+1]*Q_im;
		  tp3 += IntM[2*2]*Q_re - IntM[2*2+1]*Q_im;
		  tp4 += IntM[3*2]*Q_re - IntM[3*2+1]*Q_im;
		  tp5 += IntM[4*2]*Q_re - IntM[4*2+1]*Q_im;
		}
      }
    }	     
    //tmp0 += tp0; tmp2 += dv*tp2 + Gridv((double)j1)*tp0; tmp3 += dv*tp3 +Gridv((double)j2)*tp0 ;tmp4 += dv*tp4 + Gridv((double)j3)*tp0;  //tmp5 += (Gridv((double)j1)*Gridv((double)j1) + Gridv((double)j2)*Gridv((double)j2) + Gridv((double)j3)*Gridv((double)j3))*tp0 +dv*dv*tp5+2*dv*(Gridv((double)j1)*tp2 + Gridv((double)j2)*tp3 +Gridv((double)j3)*tp4);     
	//tmp5 += dv*dv*tp5 -  (Gridv((double)j1)*Gridv((double)j1) + Gridv((double)j2)*Gridv((double)j2) + Gridv((double)j3)*Gridv((double)j3))*tp0 + 2*(Gridv((double)j1)*(dv*tp2 + Gridv((double)j1)*tp0) + Gridv((double)j2)*(dv*tp3 + Gridv((double)j2)*tp0) + Gridv((double)j3)*(dv*tp4 + Gridv((double)j3)*tp0));
	  
    k_v = l*size_v + kt;      
    tp0 = U[k_v*6+0] + U[k_v*6+5]/4. + dt*tp0/scalev/scaleL/scale3;
    tp2 = U[k_v*6+2] + dt*tp2*12./scalev/scaleL/scale3;
    tp3 = U[k_v*6+3] + dt*tp3*12./scalev/scaleL/scale3;
    tp4 = U[k_v*6+4] + dt*tp4*12./scalev/scaleL/scale3;
    tp5 = U[k_v*6+0]/4. + U[k_v*6+5]*19./240. + dt*tp5/scalev/scaleL/scale3;

    U[k_v*6+0] = 19*tp0/4. - 15*tp5;
    U[k_v*6+5] = 60*tp5 - 15*tp0;
    U[k_v*6+2] = tp2; U[k_v*6+3] = tp3; U[k_v*6+4] = tp4;      
  }
	//printf("RK4: %g, %g, %g, %g, %g\n", tmp0, tmp2, tmp3, tmp4, 0.5*tmp5);  
}   
*/
//#else
void ComputeQ_NoMPI(double *f, fftw_complex *qHat, double **conv_weights)
{
  int i, j, k, l, m, n, x, y, z;
  int start_i, start_j, start_k, end_i, end_j, end_k;
  //fftw_complex *fftIn, *fftOut;
  double tempD, tmp0, tmp1;
  double prefactor = h_eta*h_eta*h_eta; // we dont have scale3 here.

  //fftIn = (fftw_complex *)fftw_malloc(N*N*N*sizeof(fftw_complex));
  //fftOut = (fftw_complex *)fftw_malloc(N*N*N*sizeof(fftw_complex));
  
  for(i=0;i<size_ft;i++){
    fftIn[i][0] = f[i];
    fftIn[i][1] = 0.;
  }

  fft3D(fftIn, fftOut);
  
  //printf("fft done\n");
  #pragma omp parallel for schedule(dynamic) private(j,k,l,m,n,x,y,z,start_i,start_j,start_k,end_i,end_j,end_k,tempD, tmp0, tmp1) shared(qHat, fftOut, conv_weights)
  for(i=0;i<N;i++) {
    for(j=0;j<N;j++){
      for(k=0;k<N;k++)
	{
	  //figure out the windows for the convolutions (i.e. where xi(l) and eta(i)-xi(l) are in the domain)
	  if( i < N/2 ) {
	    start_i = 0;
	    end_i = i + N/2 + 1; 
	  }
	  else {
	    start_i = i - N/2 + 1;
	    end_i = N;
	  }
	  
	  if( j < N/2 ) {
	    start_j = 0;
	    end_j = j + N/2 + 1; 
	  }
	  else {
	    start_j = j - N/2 + 1;
	    end_j = N;
	  }
	  
	  if( k < N/2 ) {
	    start_k = 0;
	    end_k = k + N/2 + 1; 
	  }
	  else {
	    start_k = k - N/2 + 1;
	    end_k = N;
	  }
	  tmp0=0.; tmp1=0.;
	  for(l=start_i;l<end_i;l++) {
	    for(m=start_j;m<end_j;m++) {
	      for(n=start_k;n<end_k;n++)
		{
		
		  x = i + N/2 - l;
		  y = j + N/2 - m;
		  z = k + N/2 - n;
		  //printf("%d, %d, %d; %d, %d, %d\n",i,j,k,l,m,n);
		 // printf("tempD=%g\n",conv_weights[k + N*(j+ N*i)][n + N*(m + N*l)]);
		  //get convolution weight		   
		  tempD = conv_weights[k + N*(j+ N*i)][n + N*(m + N*l)]; 
		  //tempD = gHat3(eta[i], eta[j], eta[k], eta[l], eta[m], eta[n]); //computing weights on the fly
		  tmp0 += prefactor*wtN[l]*wtN[m]*wtN[n]*tempD*(fftOut[n + N*(m + N*l)][0]*fftOut[z + N*(y + N*x)][0] - fftOut[n + N*(m + N*l)][1]*fftOut[z + N*(y + N*x)][1]);
		  
		  tmp1 += prefactor*wtN[l]*wtN[m]*wtN[n]*tempD*(fftOut[n + N*(m + N*l)][0]*fftOut[z + N*(y + N*x)][1] + fftOut[n + N*(m + N*l)][1]*fftOut[z + N*(y + N*x)][0]);
		 // printf(" tempD=%g, prefactor=%g, fftOut=[%g,%g] at %d, %d, %d; %d, %d, %d; %d, %d, %d\n",tempD, prefactor, fftOut[z + N*(y + N*x)][0], fftOut[z + N*(y + N*x)][1],i,j,k,l,m,n, x,y,z);
		}
	    }
	  }
	// printf("%d, %d, %d done\n", i,j,k);
	 qHat[k + N*(j + N*i)][0] = tmp0;  
	 qHat[k + N*(j + N*i)][1] = tmp1;
	 // printf("%d, %d, %d write-in done\n", i,j,k);
	}
    }
  }

 // fftw_free(fftIn);
 // fftw_free(fftOut);

}

void RK4_NoMPI(double *f, int l, fftw_complex *qHat, double **conv_weights, double *U) //4-th RK. yn=yn+(3*k1+k2+k3+k4)/6
{
  int i,j,k, j1, j2, j3, k_v, k_eta,kk;  
  double Q_re, Q_im, tp0, tp2, tp3,tp4,tp5, tmp0=0., tmp2=0., tmp3=0., tmp4=0.,tmp5=0., tem;

  FS(qHat, fftOut); 
  //ifft3D(qHat, fftOut);
  #pragma omp parallel for private(i) shared(Q,fftOut,f1,f)
  for(i=0;i<size_ft;i++){    
    Q[i] = fftOut[i][0];
    f1[i] = f[i] + dt*Q[i]*nu; //BUG: this evolution (only on node values) is not consistent with our conservation routine, which preserves the exact moments of the {1,v,|v|^{2}} approximations
  }

  ComputeQ(f1, Q1_fft, conv_weights); //collides
  conserveAllMoments(Q1_fft);   	//conserves k2	

  FS(Q1_fft, fftOut);
  //ifft3D(Q1_fft, fftOut);
  #pragma omp parallel for private(i) shared(Q,Q1,fftOut,f1,f)
  for(i=0;i<size_ft;i++){ 	
    Q1[i] = fftOut[i][0];	   
    f1[i] = f[i] +  0.5*dt*Q[i]*nu + 0.5*dt*Q1[i]*nu;
  }
  
  ComputeQ(f1, Q2_fft, conv_weights); //collides
  conserveAllMoments(Q2_fft);   //conserves k3

  FS(Q2_fft, fftOut);
  //ifft3D(Q2_fft, fftOut);
  #pragma omp parallel for private(i) shared(Q,Q1,fftOut,f1,f)
  for(i=0;i<size_ft;i++){	
    Q1[i] = fftOut[i][0];	 
    f1[i] = f[i] + 0.5*dt*Q[i]*nu + 0.5*dt*Q1[i]*nu;
  }
	 
  ComputeQ(f1, Q3_fft, conv_weights); //collides
  conserveAllMoments(Q3_fft);                //conserves k4
  
  /*FS(Q3_fft, fftOut);
  for(i=0;i<size_ft;i++){
    Q3[i] = fftOut[i][0];	  
    f[i] = f[i] + 0.5*dt*Q[i]*nu + dt*(Q1[i]+Q2[i]+Q3[i])*nu/6.;
  } */

 /* ProjectedNodeValue(qHat, Q);
  #pragma omp parallel for private(k_eta) shared(Q,f1,f)
  for(k_eta=0;k_eta<size_ft;k_eta++){  
    f1[k_eta] = f[k_eta] + dt*Q[k_eta]; //BUG: this evolution (only on node values) is not consistent with our conservation routine, which preserves the exact moments of the {1,v,|v|^{2}} approximations
  }

  ComputeQ(f1, Q1_fft, conv_weights); //collides
  conserveAllMoments(Q1_fft);   	//conserves k2	

  ProjectedNodeValue(Q1_fft, Q1);
  #pragma omp parallel for private(k_eta) shared(Q,Q1,f2,f)
  for(k_eta=0;k_eta<size_ft;k_eta++){   
    f2[k_eta] = f[k_eta] +  0.5*dt*Q[k_eta] + 0.5*dt*Q1[k_eta];
  }
  
	//printf("density(f2)=%g\n", getDensity(f2));
  ComputeQ(f2, Q2_fft, conv_weights); //collides
  conserveAllMoments(Q2_fft);   //conserves k3

  ProjectedNodeValue(Q2_fft, Q1);
  #pragma omp parallel for private(k_eta) shared(Q,Q1,f3,f)
  for(k_eta=0;k_eta<size_ft;k_eta++){		 
    f3[k_eta] = f[k_eta] + 0.5*dt*Q[k_eta] + 0.5*dt*Q1[k_eta];
  }
	 
  ComputeQ(f3, Q3_fft, conv_weights); //collides
  conserveAllMoments(Q3_fft);                //conserves k4
  */
  #pragma omp parallel for schedule(dynamic) private(j1,j2,j3,i,j,k,k_v,k_eta,kk,Q_re, Q_im, tp0, tp2,tp3,tp4, tp5) shared(l,qHat,U) //reduction(+: tmp0, tmp2, tmp3, tmp4, tmp5)
  for(int kt=0;kt<size_v;kt++){
    j3 = kt % Nv; j2 = ((kt-j3)/Nv) % Nv; j1 = (kt - j3 - Nv*j2)/(Nv*Nv);
    tp0=0.; tp2=0.; tp3=0.; tp4=0.; tp5=0.;
    for(i=0;i<N;i++){
      for(j=0;j<N;j++){
		for(k=0;k<N;k++){
		  k_eta = k + N*(j + N*i);
		  //kk = kt*size_ft+k_eta;
		  IntModes(i,j,k,j1,j2,j3,IntM); //the global IntM must be declared as threadprivate
		  
		  Q_re = nu*(0.5*qHat[k_eta][0] + (Q1_fft[k_eta][0]+Q2_fft[k_eta][0]+Q3_fft[k_eta][0])/6.);
		  Q_im = nu*(0.5*qHat[k_eta][1] + (Q1_fft[k_eta][1]+Q2_fft[k_eta][1]+Q3_fft[k_eta][1])/6.);
		  
		  //tem = scale3*wtN[i]*wtN[j]*wtN[k]*h_eta*h_eta*h_eta;
		  tp0 += IntM[0]*Q_re - IntM[1]*Q_im;
		  tp2 += IntM[1*2]*Q_re - IntM[1*2+1]*Q_im;
		  tp3 += IntM[2*2]*Q_re - IntM[2*2+1]*Q_im;
		  tp4 += IntM[3*2]*Q_re - IntM[3*2+1]*Q_im;
		  tp5 += IntM[4*2]*Q_re - IntM[4*2+1]*Q_im;
		}
      }
    }	     
    //tmp0 += tp0; tmp2 += dv*tp2 + Gridv((double)j1)*tp0; tmp3 += dv*tp3 +Gridv((double)j2)*tp0 ;tmp4 += dv*tp4 + Gridv((double)j3)*tp0;  //tmp5 += (Gridv((double)j1)*Gridv((double)j1) + Gridv((double)j2)*Gridv((double)j2) + Gridv((double)j3)*Gridv((double)j3))*tp0 +dv*dv*tp5+2*dv*(Gridv((double)j1)*tp2 + Gridv((double)j2)*tp3 +Gridv((double)j3)*tp4);     
	//tmp5 += dv*dv*tp5 -  (Gridv((double)j1)*Gridv((double)j1) + Gridv((double)j2)*Gridv((double)j2) + Gridv((double)j3)*Gridv((double)j3))*tp0 + 2*(Gridv((double)j1)*(dv*tp2 + Gridv((double)j1)*tp0) + Gridv((double)j2)*(dv*tp3 + Gridv((double)j2)*tp0) + Gridv((double)j3)*(dv*tp4 + Gridv((double)j3)*tp0));
	  
    k_v = l*size_v + kt;      
    tp0 = U[k_v*6+0] + U[k_v*6+5]/4. + dt*tp0/scalev/scaleL/scale3;
    tp2 = U[k_v*6+2] + dt*tp2*12./scalev/scaleL/scale3;
    tp3 = U[k_v*6+3] + dt*tp3*12./scalev/scaleL/scale3;
    tp4 = U[k_v*6+4] + dt*tp4*12./scalev/scaleL/scale3;
    tp5 = U[k_v*6+0]/4. + U[k_v*6+5]*19./240. + dt*tp5/scalev/scaleL/scale3;

    U[k_v*6+0] = 19*tp0/4. - 15*tp5;
    U[k_v*6+5] = 60*tp5 - 15*tp0;
    U[k_v*6+2] = tp2; U[k_v*6+3] = tp3; U[k_v*6+4] = tp4;      
  }
	//printf("RK4: %g, %g, %g, %g, %g\n", tmp0, tmp2, tmp3, tmp4, 0.5*tmp5);  
}   
//#endif

//#endif
