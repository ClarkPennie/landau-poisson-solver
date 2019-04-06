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

double S1hat_maxmols(double ki1,double ki2,double ki3)
{
	double r = sqrt(ki1*ki1+ki2*ki2+ki3*ki3);
	double Rr = R_v*r;
	double cos_Rr = cos(Rr);
	double sin_Rr = sin(Rr);
	if(r==0.) return 2*sqrt(1./(2.*PI))*pow(R_v,5.)/5.;
	else return sqrt(2./PI)*(-Rr*Rr*Rr*cos_Rr + 3*Rr*Rr*sin_Rr + 6*Rr*cos_Rr - 6*sin_Rr)/pow(r,5.);
}

double S233hat_maxmols(double ki1, double ki2, double ki3)
{
	double r = sqrt(ki1*ki1+ki2*ki2+ki3*ki3);
	double Rr = R_v*r;
	double cos_Rr = cos(Rr);
	double sin_Rr = sin(Rr);
	if(r==0.) return 2*sqrt(1./(2.*PI))*pow(R_v,5.)/15.;
	else return sqrt(2./PI)*((ki1*ki1+ki2*ki2)*(-Rr*Rr*sin_Rr - 3*Rr*cos_Rr + 3*sin_Rr) + ki3*ki3*(-Rr*Rr*Rr*cos_Rr + 5*Rr*Rr*sin_Rr + 12*Rr*cos_Rr - 12*sin_Rr))/pow(r,7.);
}

double S213hat_maxmols(double ki1, double ki2, double ki3)
{
	double r = sqrt(ki1*ki1+ki2*ki2+ki3*ki3);
	double Rr = R_v*r;
	double cos_Rr = cos(Rr);
	double sin_Rr = sin(Rr);
	if(ki1==0. || ki3==0.) return 0.;
	else return sqrt(2./PI)*ki1*ki3*(-Rr*Rr*Rr*cos_Rr + 6*Rr*Rr*sin_Rr + 15*Rr*cos_Rr - 15*sin_Rr)/(pow(r,7.));
}

double S1hat_hardspheres(double ki1,double ki2,double ki3)
{
	double r = sqrt(ki1*ki1+ki2*ki2+ki3*ki3);
	double Rr = R_v*r;
	double cos_Rr = cos(Rr);
	double sin_Rr = sin(Rr);
	if(r==0.) return sqrt(1./(2.*PI))*pow(R_v,6.)/3.;
	else return sqrt(2./PI)*(4.*(Rr*Rr - 6.)*Rr*sin_Rr - (Rr*Rr*(Rr*Rr - 12.) + 24.)*cos_Rr + 24.)/pow(r,6.);
}

double S233hat_hardspheres(double ki1, double ki2, double ki3)
{
	double r = sqrt(ki1*ki1+ki2*ki2+ki3*ki3);
	double Rr = R_v*r;
	double cos_Rr = cos(Rr);
	double sin_Rr = sin(Rr);
	if(r==0.) return sqrt(1./(2.*PI))*pow(R_v,6.)/9.;
	else return sqrt(2./PI)*((ki1*ki1+ki2*ki2)*((8. - Rr*Rr)*Rr*sin_Rr + 4.*(2. - Rr*Rr)*cos_Rr - 8.) + ki3*ki3*((Rr*Rr*(20. - Rr*Rr) - 40.)*cos_Rr + (6.*Rr*Rr - 40.)*Rr*sin_Rr + 40.))/pow(r,8.);
}

double S213hat_hardspheres(double ki1, double ki2, double ki3)
{
	double r = sqrt(ki1*ki1+ki2*ki2+ki3*ki3);
	double Rr = R_v*r;
	double cos_Rr = cos(Rr);
	double sin_Rr = sin(Rr);
	if(ki1==0. || ki3==0.) return 0.;
	else return sqrt(2./PI)*ki1*ki3*((Rr*Rr*(24. - Rr*Rr) - 48.)*cos_Rr + (7.*Rr*Rr - 48.)*Rr*sin_Rr + 48.)/(pow(r,8.));
}

double gHat3(double eta1, double eta2, double eta3, double ki1, double ki2, double ki3, int gamma)
{
	double result = 0.;
	double ki[3]={ki1,ki2,ki3}, zeta[3]={eta1,eta2,eta3}; // ki=w, zeta=xi in the notes
	double Shat[3][3];
	double r=sqrt(ki1*ki1+ki2*ki2+ki3*ki3);
	int i,j;
	
	if(gamma==-3)
	{
		Shat[0][0]=S1hat(ki1,ki2,ki3)-S233hat(ki2,ki3,ki1);
		Shat[1][1]=S1hat(ki1,ki2,ki3)-S233hat(ki1,ki3,ki2);
		Shat[2][2]=S1hat(ki1,ki2,ki3)-S233hat(ki1,ki2,ki3);
		Shat[0][1]=-S213hat(ki1,ki3,ki2);
		Shat[0][2]=-S213hat(ki1,ki2,ki3);
		Shat[1][2]=-S213hat(ki2,ki1,ki3);
		Shat[1][0]=Shat[0][1]; Shat[2][0]=Shat[0][2]; Shat[2][1]=Shat[1][2];
	}
	else if(gamma==0)
	{
		Shat[0][0]=S1hat_maxmols(ki1,ki2,ki3)-S233hat_maxmols(ki2,ki3,ki1);
		Shat[1][1]=S1hat_maxmols(ki1,ki2,ki3)-S233hat_maxmols(ki1,ki3,ki2);
		Shat[2][2]=S1hat_maxmols(ki1,ki2,ki3)-S233hat_maxmols(ki1,ki2,ki3);
		Shat[0][1]=-S213hat_maxmols(ki1,ki3,ki2);
		Shat[0][2]=-S213hat_maxmols(ki1,ki2,ki3);
		Shat[1][2]=-S213hat_maxmols(ki2,ki1,ki3);
		Shat[1][0]=Shat[0][1]; Shat[2][0]=Shat[0][2]; Shat[2][1]=Shat[1][2];
	}
	else if(gamma==1)
	{
		Shat[0][0]=S1hat_hardspheres(ki1,ki2,ki3)-S233hat_hardspheres(ki2,ki3,ki1);
		Shat[1][1]=S1hat_hardspheres(ki1,ki2,ki3)-S233hat_hardspheres(ki1,ki3,ki2);
		Shat[2][2]=S1hat_hardspheres(ki1,ki2,ki3)-S233hat_hardspheres(ki1,ki2,ki3);
		Shat[0][1]=-S213hat_hardspheres(ki1,ki3,ki2);
		Shat[0][2]=-S213hat_hardspheres(ki1,ki2,ki3);
		Shat[1][2]=-S213hat_hardspheres(ki2,ki1,ki3);
		Shat[1][0]=Shat[0][1]; Shat[2][0]=Shat[0][2]; Shat[2][1]=Shat[1][2];
	}

	if(gamma==-3)
	{
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				//result += Shat[i][j]*(2.*ki[j]-zeta[j])*zeta[i];
				result += Shat[i][j]*(zeta[i]-ki[i])*(zeta[j]-ki[j]);
			}
		}
		if(r==0.)result = -result;
		else result = (sqrt(8./PI))*(R_v*r-sin(R_v*r))/(R_v*r) - result;
	}
	else
	{
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
			  result += Shat[i][j]*(2.*ki[j]-zeta[j])*zeta[i];
			}
		}
	}
	return result;
}

double gHat_LH(double eta1_L, double eta2_L, double eta3_L, double ki1_L, double ki2_L, double ki3_L, double epsilon)
{
    double result = 0.;
    double ki_L[3]={ki1_L,ki2_L,ki3_L}, zeta_L[3]={eta1_L,eta2_L,eta3_L}; // ki=w, zeta=xi in the notes
    double Shat[3][3];
    double r=sqrt(ki1_L*ki1_L+ki2_L*ki2_L+ki3_L*ki3_L);
    double ki1_tp=epsilon*(eta1_L-ki1_L);
    double ki2_tp=epsilon*(eta2_L-ki2_L);
    double ki3_tp=epsilon*(eta3_L-ki3_L);
    
    int i,j;
    
    Shat[0][0]=S1hat(ki1_tp,ki2_tp,ki3_tp)-S233hat(ki2_tp,ki3_tp,ki1_tp);
    Shat[1][1]=S1hat(ki1_tp,ki2_tp,ki3_tp)-S233hat(ki1_tp,ki3_tp,ki2_tp);
    Shat[2][2]=S1hat(ki1_tp,ki2_tp,ki3_tp)-S233hat(ki1_tp,ki2_tp,ki3_tp);
    Shat[0][1]=-S213hat(ki1_tp,ki3_tp,ki2_tp);
    Shat[0][2]=-S213hat(ki1_tp,ki2_tp,ki3_tp);
    Shat[1][2]=-S213hat(ki2_tp,ki1_tp,ki3_tp);
    Shat[1][0]=Shat[0][1]; Shat[2][0]=Shat[0][2]; Shat[2][1]=Shat[1][2];
    
    
    for(i=0;i<3;i++){
        for(j=0;j<3;j++){
           result += Shat[i][j]*(zeta_L[i]*ki_L[j]-epsilon*epsilon*zeta_L[i]*(zeta_L[j]-ki_L[j]));  //convolution weight G in Q_LH
        }
    }
	return result;
}


double gHat_HL(double eta1_H, double eta2_H, double eta3_H, double ki1_H, double ki2_H, double ki3_H, double epsilon)
{
    double result = 0.;
    double ki_H[3]={ki1_H,ki2_H,ki3_H}, zeta_H[3]={eta1_H,eta2_H,eta3_H}; // ki=w, zeta=xi in the notes
    double Shat[3][3];
    double r=sqrt(ki1_H*ki1_H+ki2_H*ki2_H+ki3_H*ki3_H);
    double ki1_tp=ki1_H-eta1_H;
    double ki2_tp=ki2_H-eta2_H;
    double ki3_tp=ki3_H-eta3_H;
    
    int i,j;
    
    Shat[0][0]=S1hat(ki1_tp,ki2_tp,ki3_tp)-S233hat(ki2_tp,ki3_tp,ki1_tp);
    Shat[1][1]=S1hat(ki1_tp,ki2_tp,ki3_tp)-S233hat(ki1_tp,ki3_tp,ki2_tp);
    Shat[2][2]=S1hat(ki1_tp,ki2_tp,ki3_tp)-S233hat(ki1_tp,ki2_tp,ki3_tp);
    Shat[0][1]=-S213hat(ki1_tp,ki3_tp,ki2_tp);
    Shat[0][2]=-S213hat(ki1_tp,ki2_tp,ki3_tp);
    Shat[1][2]=-S213hat(ki2_tp,ki1_tp,ki3_tp);
    Shat[1][0]=Shat[0][1]; Shat[2][0]=Shat[0][2]; Shat[2][1]=Shat[1][2];

    
    for(i=0;i<3;i++){
        for(j=0;j<3;j++){
            
            result += Shat[i][j]*(zeta_H[i]*(zeta_H[j]-ki_H[j])-epsilon*epsilon*zeta_H[i]*ki_H[j]);  //convolutiong weight R in Q_HL
        }
    }
    return result;
}


double gHat3_2( double eta1, double eta2, double eta3, double ki1, double ki2, double ki3, int id)
{
	double result1 = 0., result2=0.;
	double ki[3]={ki1,ki2,ki3}, zeta[3]={eta1,eta2,eta3};// corresponding to the notes, ki=w, eta=xi
	double Shat[3][3];
	double r=sqrt(ki1*ki1+ki2*ki2+ki3*ki3);
	double eta_ki=ki1*eta1+ki2*eta2+ki3*eta3;
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
	    result2 += Shat[i][j]*zeta[j]*zeta[i];
	  }
	}
	if(r==0.)result1 = 0.;
	else result1 = sqrt(8./PI)*eta_ki*(R_v*r-sin(R_v*r))/(R_v*r*r*r);

	if(id==0) return result1;
	else return -result2;

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


void generate_conv_weights(double **conv_weights, double **conv_weights_LH, double **conv_weights_HL, int gamma)
{
  int i, j, k, l, m, n;
   #pragma omp parallel for private(i,j,k,l,m,n) shared(conv_weights)
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      for(k=0;k<N;k++){ 
	for(l=0;l<N;l++){
	  for(m=0;m<N;m++){
	    for(n=0;n<N;n++) {
	     conv_weights[k + N*(j + N*i)][n + N*(m + N*l)] = gHat3(eta[i], eta[j], eta[k], eta[l], eta[m], eta[n], gamma); // in the notes, correspondingly, (i,j,k)-kxi, (l,m,n)-w
         conv_weights_LH[k + N*(j + N*i)][n + N*(m + N*l)] = gHat_LH(eta_L[i], eta_L[j], eta_L[k], eta_L[l], eta_L[m], eta_L[n], gamma);
         conv_weights_HL[k + N*(j + N*i)][n + N*(m + N*l)] = gHat_HL(eta_H[i], eta_H[j], eta_H[k], eta_H[l], eta_H[m], eta_H[n], gamma);
	    }
	  }
	}
      }
    }
  }
}

void generate_conv_weights2(double **conv_weights, int id)
{
  int i, j, k, l, m, n;

   #pragma omp parallel for  private(i,j,k,l,m,n) shared(conv_weights)
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      for(k=0;k<N;k++){
	for(l=0;l<N;l++){
	  for(m=0;m<N;m++){
	    for(n=0;n<N;n++) {
	      conv_weights[k + N*(j + N*i)][n + N*(m + N*l)] = gHat3_2(eta[i], eta[j], eta[k], eta[l], eta[m], eta[n], id);
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
function ifft3D
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

void ComputeQ_FandL(double *f, fftw_complex *qHat, double **conv_weights, fftw_complex *qHat_linear, double **conv_weights_linear)
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

void ComputeQ(double *f, fftw_complex *qHat, double **conv_weights)
{
	int i, j, k, l, m, n, x, y, z;												// declare (i,j,k) (the indices for a given value of given ki = ki_(i,j,k)), (l,m,n) (counters for the quadrature to calculate the integral w.r.t. eta in the evaluation of qHat and so also represent the indices of a given eta = eta_(l,m,n)) & (x,y,z) (the indices for the value of a subtraction in the calculation, namely eta_(x,y,z) = ki_(i,j,k) - eta_(l,m,n))
	int start_i, start_j, start_k, end_i, end_j, end_k;							// declare start_i, start_j & start_k (the indices for the values of the lower bounds of integration in computation of the convolution, corresponding to the lowest point where both functions are non-zero, in each velocity direction) and end_i, end_j & end_k (the indices for the values of the upper bounds of integration in computation of the convolution, corresponding to the highest point where both functions are non-zero, in each velocity direction)
	double tempD, tmp0, tmp1;													// declare tempD (the value of the convolution weight at a given ki & eta), tmp0 (which will become the real part of qHat) & tmp1 (which will become the imaginary part of qHat)
	double prefactor = h_eta*h_eta*h_eta; 										// declare prefactor (the value of h_eta^3, as no scale3 in Fourier space) and set its value

    for(i=0;i<size_ft;i++)                                                        // initialise the input of the FFT
    {
        fftIn[i][0] = f[i];                                                        // set the real part to the sampling of the solution stored in f
        fftIn[i][1] = 0.;                                                        // set the imaginary part to zero
    }
    
    fft3D(fftIn, fftOut);                                                        // perform the FFT of fftIn and store the result in fftOut


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

							tempD = conv_weights[k + N*(j+ N*i)][n + N*(m + N*l)];		// set tempD to the value of the convolution weight corresponding to current value of ki(i,j,k) & eta(l,m,n)
							tmp0 += prefactor*wtN[l]*wtN[m]*wtN[n]*tempD*(fftOut[n + N*(m + N*l)][0]*fftOut[z + N*(y + N*x)][0] - fftOut[n + N*(m + N*l)][1]*fftOut[z + N*(y + N*x)][1]);		// increment the value of the real part of qHat(ki(i,j,k)) by fHat(eta(l,m,n))*fHat(ki(i,j,k)-eta(l,m,n))*conv_weight(ki(i,j,k),eta(l,m,n)) for the current values of l, m & n in the quadrature sum

							tmp1 += prefactor*wtN[l]*wtN[m]*wtN[n]*tempD*(fftOut[n + N*(m + N*l)][0]*fftOut[z + N*(y + N*x)][1] + fftOut[n + N*(m + N*l)][1]*fftOut[z + N*(y + N*x)][0]);		// increment the value of the imaginary part of qHat(ki(i,j,k)) by fHat(eta(l,m,n))*fHat(ki(i,j,k)-eta(l,m,n))*conv_weight(ki(i,j,k),eta(l,m,n)) for the current values of l, m & n in the quadrature sum

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


void ComputeQ_LH(double *f_L, double *f_H, fftw_complex *qHat, double **conv_weights_LH)
{
    int i, j, k, l, m, n, x, y, z;                                                // declare (i,j,k) (the indices for a given value of given ki = ki_(i,j,k)), (l,m,n) (counters for the quadrature to calculate the integral w.r.t. eta in the evaluation of qHat and so also represent the indices of a given eta = eta_(l,m,n)) & (x,y,z) (the indices for the value of a subtraction in the calculation, namely eta_(x,y,z) = ki_(i,j,k) - eta_(l,m,n))
    int start_i, start_j, start_k, end_i, end_j, end_k;                            // declare start_i, start_j & start_k (the indices for the values of the lower bounds of integration in computation of the convolution, corresponding to the lowest point where both functions are non-zero, in each velocity direction) and end_i, end_j & end_k (the indices for the values of the upper bounds of integration in computation of the convolution, corresponding to the highest point where both functions are non-zero, in each velocity direction)
    double tempD_LH, tmp0, tmp1;                                                    // declare tempD (the value of the convolution weight at a given ki & eta), tmp0 (which will become the real part of qHat) & tmp1 (which will become the imaginary part of qHat)
    double prefactor = h_eta_L*h_eta_L*h_eta_L;                                         // declare prefactor (the value of h_eta_L^3, as no scale3 in Fourier space) and set its value
    
    for(i=0;i<size_ft;i++)                                                        // initialise the input of the FFT
    {
        fftIn_L[i][0] = f_L[i];                                                     // set the real part to the sampling of the solution stored in f
        fftIn_L[i][1] = 0.;                                                        // set the imaginary part to zero
        fftIn_H[i][0] = f_H[i];
        fftIn_H[i][1] = 0.;
    }
    
       fft3D(fftIn_L, fftOut_L);                                                        // perform the FFT of fftIn and store the result in fftOut
       fft3D(fftIn_H, fftOut_H);
    
#pragma omp parallel for schedule(dynamic) private(i,j,k,l,m,n,x,y,z,start_i,start_j,start_k,end_i,end_j,end_k,tempD_LH) shared(qHat, fftOut_H, fftOut_L, conv_weights_LH) reduction(+:tmp0, tmp1)
    for(i=0;i<N;i++)                                                             // loop through all points in the ki_1
    {
        for(j=0;j<N;j++)                                                        // loop through all points in the ki_2
        {
            for(k=0;k<N;k++)                                                    // loop through all points in the ki_3
            {
                //figure out the windows for the convolutions (i.e. where eta(l,m,n) and ki(i,j,k)-eta(l,m,n) are in the domain)
                if( i < N/2 )                                                     // if ki_1(i) < 0 then the values of l for which f(ki_1(i)-eta_1(l))*f(eta_1(l)) give a non-zero product range need -Lv < eta_1 <= ki_1 + Lv/2, due to the support of f being -Lv to Lv
                {
                    start_i = 0;                                                // set start_i to 0 to represent -Lv as the lower bound of integration
                    end_i = i + N/2 + 1;                                        // set end_i to i + N/2 + 1 to represent eta_1(i+N/2+1/2) as the upper bound of integration (with +1 there so that all cells less than index end_i are integrated over)
                }
                else                                                             // if ki_1(i) >= 0 then the values of l for which f(ki_1(i)-eta_1(l))*f(eta_1(l)) give a non-zero product range need ki_1 - Lv/2 < eta_1 <= Lv, due to the support of f being -Lv to Lv
                {
                    start_i = i - N/2 + 1;                                        // set start_i to i - N/2 + 1 to represent eta_1(i) - Lv as the lower bound of integration (with +1 since ki_1[i-N/2] - eta_1[i] actually overlaps with eta_1[-1/2] = -Lv, where f(-Lv)=0, so start at the next index)
                    end_i = N;                                                    // set end_i to N to represent Lv as the upper bound of integration (i.e. integrate over all cells with index less than N)
                }
                
                if( j < N/2 )                                                    // if ki_2(j) < 0 then the values of m for which f(ki_2(j)-eta_2(m))*f(eta_2(m)) give a non-zero product range need -Lv < eta_2 <= ki_2 + Lv/2, due to the support of f being -Lv to Lv
                {
                    start_j = 0;                                                // set start_j to 0 to represent -Lv as the lower bound of integration
                    end_j = j + N/2 + 1;                                        // set end_j to j + N/2 + 1 to represent eta_2(j+N/2+1/2) as the upper bound of integration (with +1 there so that all cells less than index end_i are integrated over)
                }
                else                                                            // if ki_2(j) >= 0 then the values of m for which f(ki_2(j)-eta_2(m))*f(eta_2(m)) give a non-zero product range need ki_2 - Lv/2 < eta_2 <= Lv, due to the support of f being -Lv to Lv
                {
                    start_j = j - N/2 + 1;                                        // set start_j to j - N/2 + 1 to represent eta_2(j) - Lv as the lower bound of integration (with +1 since ki_2[j-N/2] - eta_2[j] actually overlaps with eta_2[-1/2] = -Lv, where f(-Lv)=0, so start at the next index)
                    end_j = N;                                                    // set end_j to N to represent Lv as the upper bound of integration (i.e. integrate over all cells with index less than N)
                }
                
                if( k < N/2 )                                                    // if ki_3(k) < 0 then the values of n for which f(ki_3(k)-eta_3(n))*f(eta_3(n)) give a non-zero product range need -Lv < eta_3 <= ki_3 + Lv/2, due to the support of f being -Lv to Lv
                {
                    start_k = 0;                                                // set start_k to 0 to represent -Lv as the lower bound of integration
                    end_k = k + N/2 + 1;                                        // set end_k to k + N/2 + 1 to represent eta_3(k+N/2+1/2) as the upper bound of integration (with +1 there so that all cells less than index end_i are integrated over)
                }
                else                                                            // if ki_3(k) >= 0 then the values of n for which f(ki_3(k)-eta_3(n))*f(eta_3(n)) give a non-zero product range need ki_3 - Lv/2 < eta_3 <= Lv, due to the support of f being -Lv to Lv
                {
                    start_k = k - N/2 + 1;                                        // set start_k to k - N/2 + 1 to represent eta_3(k) - Lv as the lower bound of integration (with +1 since ki_3[k-N/2] - eta_3[k] actually overlaps with eta_3[-1/2] = -Lv, where f(-Lv)=0, so start at the next index)
                    end_k = N;                                                    // set end_k to N to represent Lv as the upper bound of integration (i.e. integrate over all cells with index less than N)
                }
                tmp0=0.; tmp1=0.;                                                // initialise tmp0 & tmp1 at zero to begin the quadrature
                for(l=start_i;l<end_i;l++)                                        // loop through all the quadrature indices in the eta_1 direction that give non-zero contribution
                {
                    for(m=start_j;m<end_j;m++)                                    // loop through all the quadrature indices in the eta_2 direction that give non-zero contribution
                    {
                        for(n=start_k;n<end_k;n++)                                // loop through all the quadrature indices in the eta_3 direction that give non-zero contribution
                        {
                            
                            x = i + N/2 - l;                                    // set the index x to i + N/2 - l to represent the subtraction eta[x] = ki[i] - eta[l]
                            y = j + N/2 - m;                                    // set the index y to j + N/2 - m to represent the subtraction eta[y] = ki[j] - eta[m]
                            z = k + N/2 - n;                                    // set the index z to k + N/2 - n to represent the subtraction eta[z] = ki[k] - eta[n]
                            
                            tempD_LH = conv_weights_LH[k + N*(j+ N*i)][n + N*(m + N*l)];        // set tempD to the value of the convolution weight corresponding to current value of ki(i,j,k) & eta(l,m,n)
                            // MULTI-SPECIES NOTE: Will need to update fftOut terms here to _H/_L
                            tmp0 += prefactor*wtN[l]*wtN[m]*wtN[n]*tempD_LH*(fftOut_L[n + N*(m + N*l)][0]*fftOut_H[z + N*(y + N*x)][0] - fftOut_L[n + N*(m + N*l)][1]*fftOut_H[z + N*(y + N*x)][1]);        // increment the value of the real part of qHat(ki(i,j,k)) by fHat_L(eta(l,m,n))*fHat_H(eps*(ki(i,j,k)-eta(l,m,n)))*conv_weight(ki(i,j,k),eta(l,m,n)) for the current values of l, m & n in the quadrature sum
                            
                            tmp1 += prefactor*wtN[l]*wtN[m]*wtN[n]*tempD_LH*(fftOut_L[n + N*(m + N*l)][0]*fftOut_H[z + N*(y + N*x)][1] + fftOut_L[n + N*(m + N*l)][1]*fftOut_H[z + N*(y + N*x)][0]);        // increment the value of the imaginary part of qHat(ki(i,j,k)) by fHat_L(eta(l,m,n))*fHat_H(eps*(ki(i,j,k)-eta(l,m,n)))*conv_weight(ki(i,j,k),eta(l,m,n)) for the current values of l, m & n in the quadrature sum
                            
                            // printf(" tempD=%g, prefactor=%g, fftOut=[%g,%g] at %d, %d, %d; %d, %d, %d; %d, %d, %d\n",tempD, prefactor, fftOut[z + N*(y + N*x)][0], fftOut[z + N*(y + N*x)][1],i,j,k,l,m,n, x,y,z);
                        }
                    }
                }
                // printf("%d, %d, %d done\n", i,j,k);
                qHat[k + N*(j + N*i)][0] = tmp0;                                // set the real part of qHat(ki(i,j,k)) to the value tmp0 calculated in the quadrature
                qHat[k + N*(j + N*i)][1] = tmp1;                                // set the imaginary part of qHat(ki(i,j,k)) to the value tmp1 calculated in the quadrature
                // printf("%d, %d, %d write-in done\n", i,j,k);
            }
        }
    }
}


void ComputeQ_HL(double *f_L, double *f_H, fftw_complex *qHat, double **conv_weights_HL)
{
    int i, j, k, l, m, n, x, y, z;                                                // declare (i,j,k) (the indices for a given value of given ki = ki_(i,j,k)), (l,m,n) (counters for the quadrature to calculate the integral w.r.t. eta in the evaluation of qHat and so also represent the indices of a given eta = eta_(l,m,n)) & (x,y,z) (the indices for the value of a subtraction in the calculation, namely eta_(x,y,z) = ki_(i,j,k) - eta_(l,m,n))
    int start_i, start_j, start_k, end_i, end_j, end_k;                            // declare start_i, start_j & start_k (the indices for the values of the lower bounds of integration in computation of the convolution, corresponding to the lowest point where both functions are non-zero, in each velocity direction) and end_i, end_j & end_k (the indices for the values of the upper bounds of integration in computation of the convolution, corresponding to the highest point where both functions are non-zero, in each velocity direction)
    double tempD_HL, tmp0, tmp1;                                                    // declare tempD (the value of the convolution weight at a given ki & eta), tmp0 (which will become the real part of qHat) & tmp1 (which will become the imaginary part of qHat)
    double prefactor = h_eta_H*h_eta_H*h_eta_H;                                         // declare prefactor (the value of h_eta_H^3, as no scale3 in Fourier space) and set its value
    
    for(i=0;i<size_ft;i++)                                                        // initialise the input of the FFT
    {
        fftIn_L[i][0] = f_L[i];                                                        // set the real part to the sampling of the solution stored in f
        fftIn_L[i][1] = 0.;                                                        // set the imaginary part to zero
        fftIn_H[i][0] = f_H[i];
        fftIn_H[i][1] = 0.;
    }
    
    fft3D(fftIn_L, fftOut_L);                                                    // perform the FFT of fftIn and store the result in fftOut
    fft3D(fftIn_H, fftOut_H);                                                   // perform the FFT of fftIn and store the result in fftOut
    
#pragma omp parallel for schedule(dynamic) private(i,j,k,l,m,n,x,y,z,start_i,start_j,start_k,end_i,end_j,end_k,tempD_HL) shared(qHat, fftOut_L, fftOut_H, conv_weights_HL) reduction(+:tmp0, tmp1)
    for(i=0;i<N;i++)                                                             // loop through all points in the ki_1
    {
        for(j=0;j<N;j++)                                                        // loop through all points in the ki_2
        {
            for(k=0;k<N;k++)                                                    // loop through all points in the ki_3
            {
                //figure out the windows for the convolutions (i.e. where eta(l,m,n) and ki(i,j,k)-eta(l,m,n) are in the domain)
                if( i < N/2 )                                                     // if ki_1(i) < 0 then the values of l for which f(ki_1(i)-eta_1(l))*f(eta_1(l)) give a non-zero product range need -Lv < eta_1 <= ki_1 + Lv/2, due to the support of f being -Lv to Lv
                {
                    start_i = 0;                                                // set start_i to 0 to represent -Lv as the lower bound of integration
                    end_i = i + N/2 + 1;                                        // set end_i to i + N/2 + 1 to represent eta_1(i+N/2+1/2) as the upper bound of integration (with +1 there so that all cells less than index end_i are integrated over)
                }
                else                                                             // if ki_1(i) >= 0 then the values of l for which f(ki_1(i)-eta_1(l))*f(eta_1(l)) give a non-zero product range need ki_1 - Lv/2 < eta_1 <= Lv, due to the support of f being -Lv to Lv
                {
                    start_i = i - N/2 + 1;                                        // set start_i to i - N/2 + 1 to represent eta_1(i) - Lv as the lower bound of integration (with +1 since ki_1[i-N/2] - eta_1[i] actually overlaps with eta_1[-1/2] = -Lv, where f(-Lv)=0, so start at the next index)
                    end_i = N;                                                    // set end_i to N to represent Lv as the upper bound of integration (i.e. integrate over all cells with index less than N)
                }
                
                if( j < N/2 )                                                    // if ki_2(j) < 0 then the values of m for which f(ki_2(j)-eta_2(m))*f(eta_2(m)) give a non-zero product range need -Lv < eta_2 <= ki_2 + Lv/2, due to the support of f being -Lv to Lv
                {
                    start_j = 0;                                                // set start_j to 0 to represent -Lv as the lower bound of integration
                    end_j = j + N/2 + 1;                                        // set end_j to j + N/2 + 1 to represent eta_2(j+N/2+1/2) as the upper bound of integration (with +1 there so that all cells less than index end_i are integrated over)
                }
                else                                                            // if ki_2(j) >= 0 then the values of m for which f(ki_2(j)-eta_2(m))*f(eta_2(m)) give a non-zero product range need ki_2 - Lv/2 < eta_2 <= Lv, due to the support of f being -Lv to Lv
                {
                    start_j = j - N/2 + 1;                                        // set start_j to j - N/2 + 1 to represent eta_2(j) - Lv as the lower bound of integration (with +1 since ki_2[j-N/2] - eta_2[j] actually overlaps with eta_2[-1/2] = -Lv, where f(-Lv)=0, so start at the next index)
                    end_j = N;                                                    // set end_j to N to represent Lv as the upper bound of integration (i.e. integrate over all cells with index less than N)
                }
                
                if( k < N/2 )                                                    // if ki_3(k) < 0 then the values of n for which f(ki_3(k)-eta_3(n))*f(eta_3(n)) give a non-zero product range need -Lv < eta_3 <= ki_3 + Lv/2, due to the support of f being -Lv to Lv
                {
                    start_k = 0;                                                // set start_k to 0 to represent -Lv as the lower bound of integration
                    end_k = k + N/2 + 1;                                        // set end_k to k + N/2 + 1 to represent eta_3(k+N/2+1/2) as the upper bound of integration (with +1 there so that all cells less than index end_i are integrated over)
                }
                else                                                            // if ki_3(k) >= 0 then the values of n for which f(ki_3(k)-eta_3(n))*f(eta_3(n)) give a non-zero product range need ki_3 - Lv/2 < eta_3 <= Lv, due to the support of f being -Lv to Lv
                {
                    start_k = k - N/2 + 1;                                        // set start_k to k - N/2 + 1 to represent eta_3(k) - Lv as the lower bound of integration (with +1 since ki_3[k-N/2] - eta_3[k] actually overlaps with eta_3[-1/2] = -Lv, where f(-Lv)=0, so start at the next index)
                    end_k = N;                                                    // set end_k to N to represent Lv as the upper bound of integration (i.e. integrate over all cells with index less than N)
                }
                tmp0=0.; tmp1=0.;                                                // initialise tmp0 & tmp1 at zero to begin the quadrature
                for(l=start_i;l<end_i;l++)                                        // loop through all the quadrature indices in the eta_1 direction that give non-zero contribution
                {
                    for(m=start_j;m<end_j;m++)                                    // loop through all the quadrature indices in the eta_2 direction that give non-zero contribution
                    {
                        for(n=start_k;n<end_k;n++)                                // loop through all the quadrature indices in the eta_3 direction that give non-zero contribution
                        {
                            
                            x = i + N/2 - l;                                    // set the index x to i + N/2 - l to represent the subtraction eta[x] = ki[i] - eta[l]
                            y = j + N/2 - m;                                    // set the index y to j + N/2 - m to represent the subtraction eta[y] = ki[j] - eta[m]
                            z = k + N/2 - n;                                    // set the index z to k + N/2 - n to represent the subtraction eta[z] = ki[k] - eta[n]
                            
                            tempD_HL = conv_weights_HL[k + N*(j+ N*i)][n + N*(m + N*l)];        // set tempD to the value of the convolution weight corresponding to current value of ki(i,j,k) & eta(l,m,n)
                            // MULTI-SPECIES NOTE: Will need to update fftOut terms here to _H/_L
                            tmp0 += prefactor*wtN[l]*wtN[m]*wtN[n]*tempD_HL*(fftOut[n + N*(m + N*l)][0]*fftOut[z + N*(y + N*x)][0] - fftOut[n + N*(m + N*l)][1]*fftOut[z + N*(y + N*x)][1]);        // increment the value of the real part of qHat(ki(i,j,k)) by fHat(eta(l,m,n))*f(ki(i,j,k)-eta(l,m,n))*conv_weight(ki(i,j,k),eta(l,m,n)) for the current values of l, m & n in the quadrature sum
                            
                            tmp1 += prefactor*wtN[l]*wtN[m]*wtN[n]*tempD_HL*(fftOut[n + N*(m + N*l)][0]*fftOut[z + N*(y + N*x)][1] + fftOut[n + N*(m + N*l)][1]*fftOut[z + N*(y + N*x)][0]);        // increment the value of the imaginary part of qHat(ki(i,j,k)) by fHat(eta(l,m,n))*f(ki(i,j,k)-eta(l,m,n))*conv_weight(ki(i,j,k),eta(l,m,n)) for the current values of l, m & n in the quadrature sum
                            
                            // printf(" tempD=%g, prefactor=%g, fftOut=[%g,%g] at %d, %d, %d; %d, %d, %d; %d, %d, %d\n",tempD, prefactor, fftOut[z + N*(y + N*x)][0], fftOut[z + N*(y + N*x)][1],i,j,k,l,m,n, x,y,z);
                        }
                    }
                }
                // printf("%d, %d, %d done\n", i,j,k);
                qHat[k + N*(j + N*i)][0] = tmp0;                                // set the real part of qHat(ki(i,j,k)) to the value tmp0 calculated in the quadrature
                qHat[k + N*(j + N*i)][1] = tmp1;                                // set the imaginary part of qHat(ki(i,j,k)) to the value tmp1 calculated in the quadrature
                // printf("%d, %d, %d write-in done\n", i,j,k);
            }
        }
    }
}



void RK4_FandL(double *f, int l, fftw_complex *qHat, double **conv_weights, fftw_complex *qHat_linear, double **conv_weights_linear, double *U, double *dU)
{
	if(Homogeneous)
	{
		RK4_FandL_Homo(f, qHat, conv_weights, qHat_linear, conv_weights_linear, U, dU);
	}
	else
	{
		RK4_FandL_Inhomo(f, l, qHat, conv_weights, qHat_linear, conv_weights_linear, U, dU);
	}
}

void RK4(double *f, int l, fftw_complex *qHat, double **conv_weights, double *U, double *dU) //4-th RK. yn=yn+(3*k1+k2+k3+k4)/6
{
	if(Homogeneous)
	{
		RK4_Homo(f, qHat, conv_weights, U, dU);
	}
	else
	{
		RK4_Inhomo(f, l, qHat, conv_weights, U, dU);
	}
}

void RK4_FandL_Inhomo(double *f, int l, fftw_complex *qHat, double **conv_weights, fftw_complex *qHat_linear, double **conv_weights_linear, double *U, double *dU) //4-th RK. yn=yn+(3*k1+k2+k3+k4)/6
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
  conserveMoments(Q1_fft, Q1_fft_linear);   	//conserves k2
  
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
  conserveMoments(Q2_fft, Q2_fft_linear);   //conserves k3

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
  conserveMoments(Q3_fft, Q3_fft_linear);                //conserves k4

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

void RK4_Inhomo(double *f, int l, fftw_complex *qHat, double **conv_weights, double *U, double *dU) //4-th RK. yn=yn+(3*k1+k2+k3+k4)/6
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

  ComputeQ(f1, Q1_fft, conv_weights);								// calculate the Fourier transform of Q(f1,M) using conv_weights1 & conv_weights2 for the weights in the convolution, then store the results of the Fourier transform in Q1_fft
  conserveMoments(Q1_fft);   														// perform the explicit conservation calculation on Kn2^ = Q^(f1,f1) = Q1_fft

  FS(Q1_fft, fftOut);																	// set fftOut to the Fourier series representation of Q1_fft (i.e. the IFFT of Q1_fft, so that Kn^2 = fftOut = Q(Fn + dt*Kn^1, Fn + dt*Kn^1) )
  //ifft3D(Q1_fft, fftOut);
  #pragma omp parallel for private(i) shared(Q,Q1,fftOut,f1,f)
  for(i=0;i<size_ft;i++)																// calculate the second step of RK4
  {
    Q1[i] = fftOut[i][0];
    f1[i] = f[i] +  0.5*dt*Q[i]*nu + 0.5*dt*Q1[i]*nu;
  }

  ComputeQ(f1, Q2_fft, conv_weights);								// calculate the Fourier tranform of Q(f1,M) using conv_weights1 & conv_weights2 for the weights in the convolution, then store the results of the Fourier transform in Q2_fft
  conserveMoments(Q2_fft);   //conserves k3

  FS(Q2_fft, fftOut);
  //ifft3D(Q2_fft, fftOut);
  #pragma omp parallel for private(i) shared(Q,Q1,fftOut,f1,f)
  for(i=0;i<size_ft;i++)																// calculate the third step of RK4
  {
    Q1[i] = fftOut[i][0];
    f1[i] = f[i] + 0.5*Q[i]*nu + 0.5*Q1[i]*nu;
  }

  ComputeQ(f1, Q3_fft, conv_weights);								// calculate the Fourier transform of Q(f1,M) using conv_weights1 & conv_weights2 for the weights in the convolution, then store the results of the Fourier transform in Q3_fft
  conserveMoments(Q3_fft);                //conserves k4

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

void RK4_FandL_Homo(double *f, fftw_complex *qHat, double **conv_weights, fftw_complex *qHat_linear, double **conv_weights_linear, double *U, double *dU) //4-th RK. yn=yn+(3*k1+k2+k3+k4)/6
{
int i,j,k, j1, j2, j3, k_v, k_eta, kk, l_local, k_loc;
  double Q_re, Q_im, tp0, tp2, tp3,tp4,tp5, tmp0=0., tmp2=0., tmp3=0., tmp4=0.,tmp5=0., tem;

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
  conserveMoments(Q1_fft, Q1_fft_linear);   	//conserves k2

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
  conserveMoments(Q2_fft, Q2_fft_linear);   //conserves k3

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
  conserveMoments(Q3_fft, Q3_fft_linear);                //conserves k4

  #pragma omp parallel for private(i) shared(Q3_fft, Q3_fft_linear)
  for(i=0;i<size_ft;i++){
    Q3_fft[i][0] += Q3_fft_linear[i][0];
	Q3_fft[i][1] += Q3_fft_linear[i][1];
  }
	#pragma omp parallel for schedule(dynamic) private(k_loc,j1,j2,j3,i,j,k,k_v,k_eta,kk,Q_re, Q_im, tp0, tp2,tp3,tp4, tp5) shared(qHat,U, dU)   // calculate the fourth step of RK4 (still in Fourier space though?!) - reduction(+: tmp0, tmp2, tmp3, tmp4, tmp5)
	for(k_v = myrank_mpi*chunksize_dg; k_v < (myrank_mpi+1)*chunksize_dg; k_v++){
	  j3 = k_v % Nv; j2 = ((k_v-j3)/Nv) % Nv; j1 = (k_v - j3 - Nv*j2)/(Nv*Nv);
		k_loc = k_v%chunksize_dg;
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

	  tp0 = U[k_v*6+0] + U[k_v*6+5]/4. + dt*tp0/scalev/scaleL/scale3;
	  tp2 = U[k_v*6+2] + dt*tp2*12./scalev/scaleL/scale3;
	  tp3 = U[k_v*6+3] + dt*tp3*12./scalev/scaleL/scale3;
	  tp4 = U[k_v*6+4] + dt*tp4*12./scalev/scaleL/scale3;
	  tp5 = U[k_v*6+0]/4. + U[k_v*6+5]*19./240. + dt*tp5/scalev/scaleL/scale3;

	  dU[(k_loc)*5] = 19*tp0/4. - 15*tp5;
	  dU[(k_loc)*5+4] = 60*tp5 - 15*tp0;
	  dU[(k_loc)*5+1] = tp2;
	  dU[(k_loc)*5+2] = tp3;
	  dU[(k_loc)*5+3] = tp4;
	}
}

void RK4_Homo(double *f, fftw_complex *qHat, double **conv_weights, double *U, double *dU) //4-th RK. yn=yn+(3*k1+k2+k3+k4)/6
{
  int i,j,k, j1, j2, j3, k_v, k_eta, kk, k_loc;
  double Q_re, Q_im, tp0, tp2, tp3,tp4,tp5, tmp0=0., tmp2=0., tmp3=0., tmp4=0.,tmp5=0., tem;

  FS(qHat, fftOut); 																	// set fftOut to the Fourier series representation of qHat (i.e. the IFFT of qHat)
  //ifft3D(qHat, fftOut);
  #pragma omp parallel for private(i) shared(Q,fftOut,f1,f)
  for(i=0;i<size_ft;i++)																// calculate the first step of RK4
  {
    Q[i] = fftOut[i][0];																// this is Q(Fn, Fn) so that Kn^1 = dt*Q(Fn, Fn) = dt*Q[i]
    f1[i] = f[i] + dt*Q[i]*nu; 															// this is Fn + Kn^1(*nu...?) BUG: this evolution (only on node values) is not consistent with our conservation routine, which preserves the exact moments of the {1,v,|v|^{2}} approximations
  }

  ComputeQ(f1, Q1_fft, conv_weights);													// calculate the Fourier tranform of Q(f1,f1) using conv_weights for the weights in the convolution, then store the results of the Fourier transform in Q1_fft
  conserveMoments(Q1_fft);   															// perform the explicit conservation calculation on Kn2^ = Q^(f1,f1) = Q1_fft

  FS(Q1_fft, fftOut);																	// set fftOut to the Fourier series representation of Q1_fft (i.e. the IFFT of Q1_fft, so that Kn^2 = fftOut = Q(Fn + dt*Kn^1, Fn + dt*Kn^1) )
  //ifft3D(Q1_fft, fftOut);
  #pragma omp parallel for private(i) shared(Q,Q1,fftOut,f1,f)
  for(i=0;i<size_ft;i++)																// calculate the second step of RK4
  {
    Q1[i] = fftOut[i][0];
    f1[i] = f[i] +  0.5*dt*Q[i]*nu + 0.5*dt*Q1[i]*nu;
  }

  ComputeQ(f1, Q2_fft, conv_weights); //collides
  conserveMoments(Q2_fft);   //conserves k3

  FS(Q2_fft, fftOut);
  //ifft3D(Q2_fft, fftOut);
  #pragma omp parallel for private(i) shared(Q,Q1,fftOut,f1,f)
  for(i=0;i<size_ft;i++)																// calculate the third step of RK4
  {
    Q1[i] = fftOut[i][0];
    f1[i] = f[i] + 0.5*Q[i]*nu + 0.5*Q1[i]*nu;
  }

  ComputeQ(f1, Q3_fft, conv_weights);     //collides
  conserveMoments(Q3_fft);                //conserves k4

  #pragma omp parallel for schedule(dynamic) private(k_loc,j1,j2,j3,i,j,k,k_v,k_eta,kk,Q_re, Q_im, tp0, tp2,tp3,tp4, tp5) shared(qHat,U, dU)   // calculate the fourth step of RK4 (still in Fourier space though?!) - reduction(+: tmp0, tmp2, tmp3, tmp4, tmp5)
  for(k_v = myrank_mpi*chunksize_dg; k_v < (myrank_mpi+1)*chunksize_dg; k_v++){
    j3 = k_v % Nv; j2 = ((k_v-j3)/Nv) % Nv; j1 = (k_v - j3 - Nv*j2)/(Nv*Nv);
	k_loc = k_v%chunksize_dg;
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

    tp0 = U[k_v*6+0] + U[k_v*6+5]/4. + dt*tp0/scalev/scaleL/scale3;
    tp2 = U[k_v*6+2] + dt*tp2*12./scalev/scaleL/scale3;
    tp3 = U[k_v*6+3] + dt*tp3*12./scalev/scaleL/scale3;
    tp4 = U[k_v*6+4] + dt*tp4*12./scalev/scaleL/scale3;
    tp5 = U[k_v*6+0]/4. + U[k_v*6+5]*19./240. + dt*tp5/scalev/scaleL/scale3;

    dU[(k_loc)*5] = 19*tp0/4. - 15*tp5;
    dU[(k_loc)*5+4] = 60*tp5 - 15*tp0;
    dU[(k_loc)*5+1] = tp2;
    dU[(k_loc)*5+2] = tp3;
    dU[(k_loc)*5+3] = tp4;
  }
}


void RK4_Homo_L(double *f_L, double *f_H, fftw_complex *qHat_LL, fftw_complex *qHat_LH, double **conv_weights, double **conv_weights_LH, double *U, double *dU, double epsilon)
{
    int i,j,k, j1, j2, j3, k_v, k_eta, kk, k_loc;
    double Q_re, Q_im, tp0, tp2, tp3,tp4,tp5, tmp0=0., tmp2=0., tmp3=0., tmp4=0.,tmp5=0., tem;
    
    FS(qHat_LL, fftOut_LL);                                                                     // set fftOut to the Fourier series representation of qHat (i.e. the IFFT of qHat)
    FS(qHat_LH, fftOut_LH);
    
	#pragma omp parallel for private(i) shared(Q,fftOut,f1_L,f_L)
    for(i=0;i<size_ft;i++)
    {
        Q[i] = fftOut_LL[i][0] + fftOut_LH[i][0];                                                          // this is Q(Fn, Fn) so that Kn^1 = dt*Q(Fn, Fn) = dt*Q[i]
        f1_L[i] = f_L[i] + dt*Q[i]/(epsilon*epsilon);
    }
    
    ComputeQ(f1_L, Q1_fft_LL, conv_weights);
    ComputeQ_LH(f1_L, f1_H, Q1_fft_LH, conv_weights_LH);                                                   // calculate the Fourier tranform of Q(f1,f1) using conv_weights for the weights in the convolution, then store the results of the Fourier transform in Q1_fft
    
    FS(Q1_fft_LL, fftOut_LL);
    FS(Q1_fft_LH, fftOut_LH);
    
	#pragma omp parallel for private(i) shared(Q,Q1,fftOut,f1_L,f_L)
    for(i=0;i<size_ft;i++)                                                                // calculate the second step of RK4
    {
        Q1[i] = fftOut_LL[i][0] + fftOut_LH[i][0];
        f1_L[i] = f_L[i] + 0.5*dt*Q[i]/(epsilon*epsilon) + 0.5*dt*Q1[i]/(epsilon*epsilon);
    }
    
    ComputeQ(f1_L, Q2_fft_LL, conv_weights);
    ComputeQ_LH(f1_L, f1_H, Q2_fft_LH, conv_weights_LH);
   
    FS(Q2_fft_LL, fftOut_LL);
    FS(Q2_fft_LH, fftOut_LH);
   
	#pragma omp parallel for private(i) shared(Q,Q1,fftOut,f1_L,f_L)
    for(i=0;i<size_ft;i++)                                                                // calculate the third step of RK4
    {
        Q1[i] = fftOut_LL[i][0] + fftOut_LH[i][0];
        f1_L[i] = f_L[i] + 0.5*Q[i]/(epsilon*epsilon) + 0.5*Q1[i]/(epsilon*epsilon);
    }
    
    ComputeQ(f1_L, Q3_fft_LL, conv_weights);
    ComputeQ_LH(f1_L, f1_H, Q3_fft_LH, conv_weights_LH);

    // MULTI-SPECIES NOTE: Will probably need to change which U coefficients are being used here (most likely U_L?)
	#pragma omp parallel for schedule(dynamic) private(k_loc,j1,j2,j3,i,j,k,k_v,k_eta,kk,Q_re, Q_im, tp0, tp2,tp3,tp4, tp5) shared(qHat_LL, qHat_LH,U, dU)   // calculate the fourth step of RK4 (still in Fourier space though?!) - reduction(+: tmp0, tmp2, tmp3, tmp4, tmp5)
    for(k_v = myrank_mpi*chunksize_dg; k_v < (myrank_mpi+1)*chunksize_dg; k_v++){
        j3 = k_v % Nv; j2 = ((k_v-j3)/Nv) % Nv; j1 = (k_v - j3 - Nv*j2)/(Nv*Nv);
        k_loc = k_v%chunksize_dg;
        tp0=0.; tp2=0.; tp3=0.; tp4=0.; tp5=0.;
        for(i=0;i<N;i++){
            for(j=0;j<N;j++){
                for(k=0;k<N;k++){
                    k_eta = k + N*(j + N*i);
                    
                    IntModes(i,j,k,j1,j2,j3,IntM); //the global IntM must be declared as threadprivate; BUG: forgot to uncomment this and thus IntM=0 !!
                    
                    Q_re = (0.5*qHat_LL[k_eta][0] + (Q1_fft_LL[k_eta][0]+Q2_fft_LL[k_eta][0]+Q3_fft_LL[k_eta][0])/6.)/(epsilon*epsilon) + (0.5*qHat_LH[k_eta][0] + (Q1_fft_LH[k_eta][0]+Q2_fft_LH[k_eta][0]+Q3_fft_LH[k_eta][0])/6.)/(epsilon*epsilon);
                     Q_re = (0.5*qHat_LL[k_eta][1] + (Q1_fft_LL[k_eta][1]+Q2_fft_LL[k_eta][1]+Q3_fft_LL[k_eta][1])/6.)/(epsilon*epsilon) + (0.5*qHat_LH[k_eta][1] + (Q1_fft_LH[k_eta][1]+Q2_fft_LH[k_eta][1]+Q3_fft_LH[k_eta][1])/6.)/(epsilon*epsilon);
                    
                    //tem = scale3*wtN[i]*wtN[j]*wtN[k]*h_eta*h_eta*h_eta;
                    tp0 += IntM[0]*Q_re - IntM[1]*Q_im;
                    tp2 += IntM[1*2]*Q_re - IntM[1*2+1]*Q_im;
                    tp3 += IntM[2*2]*Q_re - IntM[2*2+1]*Q_im;
                    tp4 += IntM[3*2]*Q_re - IntM[3*2+1]*Q_im;
                    tp5 += IntM[4*2]*Q_re - IntM[4*2+1]*Q_im;
                }
            }
        }
        
        tp0 = U[k_v*6+0] + U[k_v*6+5]/4. + dt*tp0/scalev/scaleL/scale3;
        tp2 = U[k_v*6+2] + dt*tp2*12./scalev/scaleL/scale3;
        tp3 = U[k_v*6+3] + dt*tp3*12./scalev/scaleL/scale3;
        tp4 = U[k_v*6+4] + dt*tp4*12./scalev/scaleL/scale3;
        tp5 = U[k_v*6+0]/4. + U[k_v*6+5]*19./240. + dt*tp5/scalev/scaleL/scale3;
        
        dU[(k_loc)*5] = 19*tp0/4. - 15*tp5;
        dU[(k_loc)*5+4] = 60*tp5 - 15*tp0;
        dU[(k_loc)*5+1] = tp2;
        dU[(k_loc)*5+2] = tp3;
        dU[(k_loc)*5+3] = tp4;
    }
}


void RK4_Homo_H(double *f_L, double *f_H, fftw_complex *qHat_HH, fftw_complex *qHat_HL, double **conv_weights, double **conv_weights_HL, double *U, double *dU, double epsilon)
{
    int i,j,k, j1, j2, j3, k_v, k_eta, kk, k_loc;
    double Q_re, Q_im, tp0, tp2, tp3,tp4,tp5, tmp0=0., tmp2=0., tmp3=0., tmp4=0.,tmp5=0., tem;
    
    FS(qHat_HH, fftOut_HH);                                                                     // set fftOut to the Fourier series representation of qHat (i.e. the IFFT of qHat)
    FS(qHat_HL, fftOut_HL);
    
	#pragma omp parallel for private(i) shared(Q,fftOut,f1_H,f_H)
    for(i=0;i<size_ft;i++)
    {
        Q[i] = fftOut_HH[i][0] + fftOut_HL[i][0];                                                          // this is Q(Fn, Fn) so that Kn^1 = dt*Q(Fn, Fn) = dt*Q[i]
        f1_H[i] = f_H[i] + dt*Q[i]/epsilon;
    }
    
    ComputeQ(f1_H, Q1_fft_HH, conv_weights);
    ComputeQ_HL(f1_L, f1_H, Q1_fft_HL, conv_weights_HL);                                                   // calculate the Fourier tranform of Q(f1,f1) using conv_weights for the weights in the convolution, then store the results of the Fourier transform in Q1_fft
    
    FS(Q1_fft_HH, fftOut_HH);
    FS(Q1_fft_HL, fftOut_HL);
    
	#pragma omp parallel for private(i) shared(Q,Q1,fftOut,f1_H,f_H)
    for(i=0;i<size_ft;i++)                                                                // calculate the second step of RK4
    {
        Q1[i] = fftOut_HH[i][0] + fftOut_HL[i][0];
        f1_H[i] = f_H[i] + 0.5*dt*Q[i]/epsilon + 0.5*dt*Q1[i]/epsilon;
    }
    
    ComputeQ(f1_H, Q2_fft_HH, conv_weights);
    ComputeQ_HL(f1_L, f1_H, Q2_fft_HL, conv_weights_HL);
    
    FS(Q2_fft_HH, fftOut_HH);
    FS(Q2_fft_HL, fftOut_HL);
    
	#pragma omp parallel for private(i) shared(Q,Q1,fftOut,f1_H,f_H)
    for(i=0;i<size_ft;i++)                                                                // calculate the third step of RK4
    {
        Q1[i] = fftOut_HH[i][0] + fftOut_HL[i][0];
        f1_H[i] = f_H[i] + 0.5*Q[i]/epsilon + 0.5*Q1[i]/epsilon;
    }
    
    ComputeQ(f1_H, Q3_fft_HH, conv_weights);
    ComputeQ_HL(f1_L, f1_H, Q3_fft_HL, conv_weights_HL);
    
    
	#pragma omp parallel for schedule(dynamic) private(k_loc,j1,j2,j3,i,j,k,k_v,k_eta,kk,Q_re, Q_im, tp0, tp2,tp3,tp4, tp5) shared(qHat_HH,qHat_HL,U, dU)   // calculate the fourth step of RK4 (still in Fourier space though?!) - reduction(+: tmp0, tmp2, tmp3, tmp4, tmp5)
    for(k_v = myrank_mpi*chunksize_dg; k_v < (myrank_mpi+1)*chunksize_dg; k_v++){
        j3 = k_v % Nv; j2 = ((k_v-j3)/Nv) % Nv; j1 = (k_v - j3 - Nv*j2)/(Nv*Nv);
        k_loc = k_v%chunksize_dg;
        tp0=0.; tp2=0.; tp3=0.; tp4=0.; tp5=0.;
        for(i=0;i<N;i++){
            for(j=0;j<N;j++){
                for(k=0;k<N;k++){
                    k_eta = k + N*(j + N*i);
                    
                    IntModes(i,j,k,j1,j2,j3,IntM); //the global IntM must be declared as threadprivate; BUG: forgot to uncomment this and thus IntM=0 !!
                    
                    Q_re = (0.5*qHat_HH[k_eta][0] + (Q1_fft_HH[k_eta][0]+Q2_fft_HH[k_eta][0]+Q3_fft_HH[k_eta][0])/6.)/epsilon + (0.5*qHat_HL[k_eta][0] + (Q1_fft_HL[k_eta][0]+Q2_fft_HL[k_eta][0]+Q3_fft_HL[k_eta][0])/6.)/epsilon;
                    Q_im = (0.5*qHat_HH[k_eta][1] + (Q1_fft_HH[k_eta][1]+Q2_fft_HH[k_eta][1]+Q3_fft_HH[k_eta][1])/6.)/epsilon + (0.5*qHat_HL[k_eta][1] + (Q1_fft_HL[k_eta][1]+Q2_fft_HL[k_eta][1]+Q3_fft_HL[k_eta][1])/6.)/epsilon;
                    
                    //tem = scale3*wtN[i]*wtN[j]*wtN[k]*h_eta*h_eta*h_eta;
                    tp0 += IntM[0]*Q_re - IntM[1]*Q_im;
                    tp2 += IntM[1*2]*Q_re - IntM[1*2+1]*Q_im;
                    tp3 += IntM[2*2]*Q_re - IntM[2*2+1]*Q_im;
                    tp4 += IntM[3*2]*Q_re - IntM[3*2+1]*Q_im;
                    tp5 += IntM[4*2]*Q_re - IntM[4*2+1]*Q_im;
                }
            }
        }
        
        tp0 = U[k_v*6+0] + U[k_v*6+5]/4. + dt*tp0/scalev/scaleL/scale3;
        tp2 = U[k_v*6+2] + dt*tp2*12./scalev/scaleL/scale3;
        tp3 = U[k_v*6+3] + dt*tp3*12./scalev/scaleL/scale3;
        tp4 = U[k_v*6+4] + dt*tp4*12./scalev/scaleL/scale3;
        tp5 = U[k_v*6+0]/4. + U[k_v*6+5]*19./240. + dt*tp5/scalev/scaleL/scale3;
        
        dU[(k_loc)*5] = 19*tp0/4. - 15*tp5;
        dU[(k_loc)*5+4] = 60*tp5 - 15*tp0;
        dU[(k_loc)*5+1] = tp2;
        dU[(k_loc)*5+2] = tp3;
        dU[(k_loc)*5+3] = tp4;
    }
}


void ComputeDFTofMaxwellian(double *UMaxwell, double **fMaxwell, fftw_complex **DFTMax)				// function to compute the Fourier transform of the initial Maxwellian
{
	setInit_spectral(UMaxwell, fMaxwell); 												// Take the coefficient of the DG solution from the advection step, and project them onto the grid used for the spectral method to perform the collision step

	for(int l=chunk_Nx*myrank_mpi;l<chunk_Nx*(myrank_mpi+1) && l<Nx;l++)
	{
		for(int i=0;i<size_ft;i++)														// initialise the input of the FFT
		{
			fftIn[i][0] = fMaxwell[l%chunk_Nx][i];												// set the real part to the sampling of the Maxwellian stored in fMaxwell
			fftIn[i][1] = 0.;														// set the imaginary part to zero
			fft3D(fftIn, DFTMax[l%chunk_Nx]);														// perform the FFT of fftIn and store the result in DFTMaxwell
		}
	}

}

void ComputeQLinear(double *f, fftw_complex *Maxwell_fftOut, fftw_complex *qHat, double **conv_weights)		// function to calculate the discrete Fourier transform of the linear collision operator Q(f, M) at a given f(Hat)
{
	int i, j, k, l, m, n, x, y, z;												// declare (i,j,k) (the indices for a given value of given ki = ki_(i,j,k)), (l,m,n) (counters for the quadrature to calculate the integral w.r.t. eta in the evaluation of qHat and so also represent the indices of a given eta = eta_(l,m,n)) & (x,y,z) (the indices for the value of a subtraction in the calculation, namely eta_(x,y,z) = ki_(i,j,k) - eta_(l,m,n))
	int start_i, start_j, start_k, end_i, end_j, end_k;							// declare start_i, start_j & start_k (the indices for the values of the lower bounds of integration in computation of the convolution, corresponding to the lowest point where both functions are non-zero, in each velocity direction) and end_i, end_j & end_k (the indices for the values of the upper bounds of integration in computation of the convolution, corresponding to the highest point where both functions are non-zero, in each velocity direction)
	double tempD, tmp0, tmp1;													// declare tempD (the value of the convolution weight at a given ki & eta), tmp0 (which will become the real part of qHat) & tmp1 (which will become the imaginary part of qHat)
//	double tempD1, tempD2;														// declare tempD1 & tempD2 (the values of the convolution weights at a given ki & eta, separated to allow for multispecies calculations with different masses)
	double prefactor = h_eta*h_eta*h_eta; 										// declare prefactor (the value of h_eta^3, as no scale3 in Fourier space) and set its value

	for(i=0;i<size_ft;i++)														// initialise the input of the FFT
	{
		fftIn[i][0] = f[i];														// set the real part to the sampling of the solution stored in f
		fftIn[i][1] = 0.;														// set the imaginary part to zero
	}

	fft3D(fftIn, fftOut);														// perform the FFT of fftIn and store the result in fftOut

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

							tempD = conv_weights[k + N*(j+ N*i)][n + N*(m + N*l)];
//							tempD1 = conv_weightsA[k + N*(j+ N*i)][n + N*(m + N*l)];
//							tempD2 = conv_weightsB[k + N*(j+ N*i)][n + N*(m + N*l)];
							//multiply the weighted fourier coeff product
							tmp0 += prefactor*wtN[l]*wtN[m]*wtN[n]*(tempD*(Maxwell_fftOut[n + N*(m + N*l)][0]*fftOut[z + N*(y + N*x)][0] - Maxwell_fftOut[n + N*(m + N*l)][1]*fftOut[z + N*(y + N*x)][1]));
							tmp1 += prefactor*wtN[l]*wtN[m]*wtN[n]*(tempD*(Maxwell_fftOut[n + N*(m + N*l)][0]*fftOut[z + N*(y + N*x)][1] + Maxwell_fftOut[n + N*(m + N*l)][1]*fftOut[z + N*(y + N*x)][0]));
//							tmp0 += prefactor*wtN[l]*wtN[m]*wtN[n]*((2*tempD1+tempD2)*(Maxwell_fftOut[n + N*(m + N*l)][0]*fftOut[z + N*(y + N*x)][0] - Maxwell_fftOut[n + N*(m + N*l)][1]*fftOut[z + N*(y + N*x)][1]));
//							tmp1 += prefactor*wtN[l]*wtN[m]*wtN[n]*((2*tempD1+tempD2)*(Maxwell_fftOut[n + N*(m + N*l)][0]*fftOut[z + N*(y + N*x)][1] + Maxwell_fftOut[n + N*(m + N*l)][1]*fftOut[z + N*(y + N*x)][0]));
						}
					}
				}
				qHat[k + N*(j + N*i)][0] = tmp0;								// set the real part of qHat(ki(i,j,k)) to the value tmp0 calculated in the quadrature
				qHat[k + N*(j + N*i)][1] = tmp1;								// set the imaginary part of qHat(ki(i,j,k)) to the value tmp1 calculated in the quadrature
			}
		}
	}
}

//void RK4Linear(double *f, fftw_complex *MaxwellHat, int l, double nu_val, fftw_complex *qHat, double **conv_weights, double *U, double *dU) //4-th RK. yn=yn+(3*k1+k2+k3+k4)/6 for the linear collision operator with a Maxwellian
void RK4Linear(double *f, fftw_complex *MaxwellHat, int l, fftw_complex *qHat, double **conv_weights, double *U, double *dU)
{
  double nu_val = nu;	// temp without nu vector
  int i,j,k, j1, j2, j3, k_v, k_eta, kk,l_local;
  double Q_re, Q_im, tp0, tp2, tp3,tp4,tp5, tmp0=0., tmp2=0., tmp3=0., tmp4=0.,tmp5=0., tem;

  l_local = l%chunk_Nx;

  FS(qHat, fftOut); 																	// set fftOut to the Fourier series representation of qHat (i.e. the IFFT of qHat)
  //ifft3D(qHat, fftOut);
  #pragma omp parallel for private(i) shared(Q,fftOut,f1,f)
  for(i=0;i<size_ft;i++)																// calculate the first step of RK4
  {
    Q[i] = fftOut[i][0];																// this is Q(Fn, Fn) so that Kn^1 = dt*Q(Fn, Fn) = dt*Q[i]
    f1[i] = f[i] + dt*Q[i]*nu_val; 															// this is Fn + Kn^1(*nu...?) BUG: this evolution (only on node values) is not consistent with our conservation routine, which preserves the exact moments of the {1,v,|v|^{2}} approximations
  }

  ComputeQLinear(f1, MaxwellHat, Q1_fft, conv_weights);								// calculate the Fourier transform of Q(f1,M) using conv_weights1 & conv_weights2 for the weights in the convolution, then store the results of the Fourier transform in Q1_fft
  conserveMoments(Q1_fft);   														// perform the explicit conservation calculation on Kn2^ = Q^(f1,f1) = Q1_fft

  FS(Q1_fft, fftOut);																	// set fftOut to the Fourier series representation of Q1_fft (i.e. the IFFT of Q1_fft, so that Kn^2 = fftOut = Q(Fn + dt*Kn^1, Fn + dt*Kn^1) )
  //ifft3D(Q1_fft, fftOut);
  #pragma omp parallel for private(i) shared(Q,Q1,fftOut,f1,f)
  for(i=0;i<size_ft;i++)																// calculate the second step of RK4
  {
    Q1[i] = fftOut[i][0];
    f1[i] = f[i] +  0.5*dt*Q[i]*nu_val + 0.5*dt*Q1[i]*nu_val;
  }

  ComputeQLinear(f1, MaxwellHat, Q2_fft, conv_weights);								// calculate the Fourier tranform of Q(f1,M) using conv_weights1 & conv_weights2 for the weights in the convolution, then store the results of the Fourier transform in Q2_fft
  conserveMoments(Q2_fft);   //conserves k3

  FS(Q2_fft, fftOut);
  //ifft3D(Q2_fft, fftOut);
  #pragma omp parallel for private(i) shared(Q,Q1,fftOut,f1,f)
  for(i=0;i<size_ft;i++)																// calculate the third step of RK4
  {
    Q1[i] = fftOut[i][0];
    f1[i] = f[i] + 0.5*Q[i]*nu_val + 0.5*Q1[i]*nu_val;
  }

  ComputeQLinear(f1, MaxwellHat, Q3_fft, conv_weights);								// calculate the Fourier transform of Q(f1,M) using conv_weights1 & conv_weights2 for the weights in the convolution, then store the results of the Fourier transform in Q3_fft
  conserveMoments(Q3_fft);                //conserves k4

  #pragma omp parallel for schedule(dynamic) private(j1,j2,j3,i,j,k,k_v,k_eta,kk,Q_re, Q_im) shared(l, l_local, qHat,U, dU) reduction(+:tp0, tp2,tp3,tp4, tp5)  // calculate the fourth step of RK4 (still in Fourier space though?!) - reduction(+: tmp0, tmp2, tmp3, tmp4, tmp5)
  for(int kt=0;kt<size_v;kt++){
    j3 = kt % Nv; j2 = ((kt-j3)/Nv) % Nv; j1 = (kt - j3 - Nv*j2)/(Nv*Nv);
    tp0=0.; tp2=0.; tp3=0.; tp4=0.; tp5=0.;
    for(i=0;i<N;i++){
      for(j=0;j<N;j++){
		for(k=0;k<N;k++){
		  k_eta = k + N*(j + N*i);

		  IntModes(i,j,k,j1,j2,j3,IntM); //the global IntM must be declared as threadprivate; BUG: forgot to uncomment this and thus IntM=0 !!

		  Q_re = nu_val*(0.5*qHat[k_eta][0] + (Q1_fft[k_eta][0]+Q2_fft[k_eta][0]+Q3_fft[k_eta][0])/6.);
		  Q_im = nu_val*(0.5*qHat[k_eta][1] + (Q1_fft[k_eta][1]+Q2_fft[k_eta][1]+Q3_fft[k_eta][1])/6.);

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

