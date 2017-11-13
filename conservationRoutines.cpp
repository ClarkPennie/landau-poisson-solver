//#include "conservationRoutines.h"
extern fftw_plan p_forward; 
extern fftw_plan p_backward; 
extern fftw_complex *temp;

void solveWithCCt(int nElem, double *b)
{
	double tmp, sum = 0.0, btmp[nElem];
	int i, j, k;
	
	/*for(k=0;k<nElem-1;k++)
	{
		i = pivot[k];
		if(pivot[k] != k)
		{
			temp = b[i];
			b[i] = b[k];
			b[k] = temp;
		}
		
		for(i=k+1;i<nElem;i++)
		b[i] -= CCt[i][k]*b[k];
	}
	
	b[nElem-1] = b[nElem-1]/CCt[nElem-1][nElem-1];
	
	for(i=nElem-2;i>=0;i--)
	{
		sum = 0.0;
		for(j=i+1;j<nElem;j++)
		sum += CCt[i][j]*b[j];
		
		b[i] = 1.0/CCt[i][i]*(b[i] - sum);
	} */
	for(i=0;i<nElem;i++)btmp[i]=b[i];
	
	if(nElem==2){ //specially when only have conservation of mass and energy
	  for(i=0;i<nElem;i++){
	    tmp=0.;
	    for(j=0;j<nElem;j++){
		    tmp += CCt_linear[j+i*nElem]*btmp[j];
	    }
	    b[i] = tmp;		
	  }	    
	}
	else
	{	  
	  for(i=0;i<nElem;i++){
	    tmp=0.;
	    for(j=0;j<nElem;j++){
		    tmp += CCt[j+i*nElem]*btmp[j];
	    }
	    b[i] = tmp;		
	  }
	}
}	

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

//MAIN CONSERVATION ROUTINE

#ifdef FullandLinear
void conserveAllMoments(fftw_complex *qHat, fftw_complex *qHat_linear) // Q^ -->Q--->conserve-->new Q^
{
	//double  *lamb, lamb_linear[2];
	double tmp0=0., tmp1=0., tmp2=0., tmp3=0., tmp4=0., tp0=0., tp1=0.;
	int i,k_eta;	
	
	//lamb = (double *)malloc(M*sizeof(double));
	//for(i=0;i<M;i++) lamb[i]=0.;
	
	#pragma omp parallel for private(k_eta) shared(qHat, qHat_linear) reduction(+:tmp0,tmp1,tmp2,tmp3,tmp4, tp0, tp1)
	for(k_eta=0;k_eta<size_ft;k_eta++){	  
	  tmp0 += qHat[k_eta][0]*C1[0][k_eta];
	  tmp1 += qHat[k_eta][1]*C2[1][k_eta];
	  tmp2 += qHat[k_eta][1]*C2[2][k_eta];
	  tmp3 += qHat[k_eta][1]*C2[3][k_eta];
	  tmp4 += qHat[k_eta][0]*C1[4][k_eta]; 

      tp0 += qHat_linear[k_eta][0]*C1[0][k_eta];	tp1 += qHat_linear[k_eta][0]*C1[4][k_eta];	  
	}
	lamb[0]=tmp0; lamb[1]=tmp1; lamb[2]=tmp2; lamb[3]=tmp3;  lamb[4]=tmp4; 	
	solveWithCCt(M, lamb); 
	
	lamb_linear[0]=tp0; lamb_linear[1]=tp1;
	solveWithCCt(2, lamb_linear); 
	//tmp0=0.; tmp1=0.;
	#pragma omp parallel for private(k_eta) shared(lamb, qHat, qHat_linear) //reduction(+:tmp0,tmp1)
	for(k_eta=0;k_eta<size_ft;k_eta++){	  
	  qHat[k_eta][0] -= (C1[0][k_eta]*lamb[0] + C1[4][k_eta]*lamb[4]);
	  qHat[k_eta][1] -= (C2[1][k_eta]*lamb[1] + C2[2][k_eta]*lamb[2] + C2[3][k_eta]*lamb[3]);	

	  qHat_linear[k_eta][0] -= (C1[0][k_eta]*lamb_linear[0] + C1[4][k_eta]*lamb_linear[1]);
	}	

	//free(lamb); 	
}
#else
void conserveAllMoments(fftw_complex *qHat) // Q^ -->Q--->conserve-->new Q^
{
	//double  *lamb, lamb_linear[2];
	double tmp0=0., tmp1=0., tmp2=0., tmp3=0., tmp4=0., tp0=0., tp1=0.;
	int i,k_eta;	
	
	//lamb = (double *)malloc(M*sizeof(double));
	//for(i=0;i<M;i++) lamb[i]=0.;
	
	#pragma omp parallel for private(k_eta) shared(qHat) reduction(+:tmp0,tmp1,tmp2,tmp3,tmp4, tp0, tp1)
	for(k_eta=0;k_eta<size_ft;k_eta++)							// calculate C*Q^
	{
	  tmp0 += qHat[k_eta][0]*C1[0][k_eta];
	  tmp1 += qHat[k_eta][1]*C2[1][k_eta];
	  tmp2 += qHat[k_eta][1]*C2[2][k_eta];
	  tmp3 += qHat[k_eta][1]*C2[3][k_eta];
	  tmp4 += qHat[k_eta][0]*C1[4][k_eta];     
	}
	lamb[0]=tmp0; lamb[1]=tmp1; lamb[2]=tmp2; lamb[3]=tmp3;  lamb[4]=tmp4; 	
	solveWithCCt(M, lamb); 										// calculate (C*C^T)^(-1)*lamb = (C*C^T)^(-1)*C*Q^ (i.e. solve (C*C^T)*X=lamb for X)
	
	#pragma omp parallel for private(k_eta) shared(lamb, qHat) //reduction(+:tmp0,tmp1)
	for(k_eta=0;k_eta<size_ft;k_eta++)							// calculate Q^ - (C^T)*(C*C^T)^(-1)*C*Q^ (the conserved Q^)
	{
	  qHat[k_eta][0] -= (C1[0][k_eta]*lamb[0] + C1[4][k_eta]*lamb[4]);
	  qHat[k_eta][1] -= (C2[1][k_eta]*lamb[1] + C2[2][k_eta]*lamb[2] + C2[3][k_eta]*lamb[3]);		 
	}	

	//free(lamb); 	
}
#endif
/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

void createCCtAndPivot()
{
	double det = 1, tmp, tmp1, tmp2;
	double a1, a2, a3;
	int i, j, k, k_eta, m1,m2,m3, iError = 0;
	int i1, i2, i3, j1, j2, j3;
	
	int pivotArray[M];
	double lapackWorkspace[M*M];
	int errinfo,nele=M,lwork=M*M;
	
	int pivotArray1[2];
	double lapackWorkspace1[4];
	int nele1=2,lwork1=4;

	#pragma omp parallel for private(a1,a2,a3,k_eta,i,j,k) shared(C1, C2)
	for(k_eta=0;k_eta<size_ft;k_eta++){	
	  k = k_eta%N;  j = ((k_eta-k)%(N*N))/N;  i = (k_eta-k-j*N)/(N*N);
	  //tmp = wtN[i]*wtN[j]*wtN[k];
	  C1[0][k_eta] = sinc(L_v*eta[i])*sinc(L_v*eta[j])*sinc(L_v*eta[k]);	
	  C1[1][k_eta] = 0.;	
	  C1[2][k_eta] = 0.;	
   	  C1[3][k_eta] = 0.;	  

	  if(eta[i] != 0) a1 = ((eta[i]*eta[i]*L_v*L_v-2)*sin(eta[i]*L_v) + 2*eta[i]*L_v*cos(eta[i]*L_v))/(eta[i]*eta[i]*eta[i]*L_v);//  sinc(eta[i]*L_v)*L_v*L_v - 2*(sinc(eta[i]*L_v) - cos(eta[i]*L_v))/eta[i]/eta[i];
	  else a1 = L_v*L_v/3.;
	  if(eta[j] != 0) a2 = ((eta[j]*eta[j]*L_v*L_v-2)*sin(eta[j]*L_v) + 2*eta[j]*L_v*cos(eta[j]*L_v))/(eta[j]*eta[j]*eta[j]*L_v);//sinc(eta[j]*L_v)*L_v*L_v - 2*(sinc(eta[j]*L_v) - cos(eta[j]*L_v))/eta[j]/eta[j];
	  else a2 = L_v*L_v/3.;
	  if(eta[k] != 0) a3 = ((eta[k]*eta[k]*L_v*L_v-2)*sin(eta[k]*L_v) + 2*eta[k]*L_v*cos(eta[k]*L_v))/(eta[k]*eta[k]*eta[k]*L_v);//sinc(eta[k]*L_v)*L_v*L_v - 2*(sinc(eta[k]*L_v) - cos(eta[k]*L_v))/eta[k]/eta[k];
	  else a3 = L_v*L_v/3.;
	 
	  C1[4][k_eta] = 0.5*(a1*sinc(L_v*eta[j])*sinc(L_v*eta[k]) + a2*sinc(L_v*eta[i])*sinc(L_v*eta[k]) + a3*sinc(L_v*eta[i])*sinc(L_v*eta[j]));
	  
	  C2[0][k_eta] = 0.;	  
	  if(eta[i]!=0) C2[1][k_eta] = -sinc(L_v*eta[j])*sinc(L_v*eta[k])*(sinc(eta[i]*L_v) - cos(eta[i]*L_v))/eta[i];	
	  else C2[1][k_eta] = 0.;	  
	  if(eta[j]!=0) C2[2][k_eta] = -sinc(L_v*eta[i])*sinc(L_v*eta[k])*(sinc(eta[j]*L_v) - cos(eta[j]*L_v))/eta[j];	
	  else C2[2][k_eta] = 0.;	  
	  if(eta[k]!=0) C2[3][k_eta] = -sinc(L_v*eta[i])*sinc(L_v*eta[j])*(sinc(eta[k]*L_v) - cos(eta[k]*L_v))/eta[k];
	  else C2[3][k_eta] = 0.;	  
	  C2[4][k_eta] = 0.; 	
	} 
	
	for(i=0;i<M;i++){
	for(j=0;j<M;j++){
		tmp = 0.;
		for(k=0;k<size_ft;k++) tmp += C1[i][k]*C1[j][k] + C2[i][k]*C2[j][k];
		CCt[i*M+j] = tmp;
		//printf("%g ", CCt[i*M+j]);
	}
	//printf("\n");
	}

	dgetrf(&nele,&nele,CCt,&nele,pivotArray,&errinfo); 
	dgetri(&nele,CCt,&nele,pivotArray,lapackWorkspace,&lwork,&errinfo);
	
	#ifdef FullandLinear
	tmp = 0.; tmp1=0.; tmp2=0.; 
	for(k=0;k<size_ft;k++) {
	   tmp += C1[0][k]*C1[0][k] + C2[0][k]*C2[0][k]; 
	   tmp1 += C1[0][k]*C1[4][k] + C2[0][k]*C2[4][k];
	   tmp2 += C1[4][k]*C1[4][k] + C2[4][k]*C2[4][k];
	}
	CCt_linear[0] = tmp; CCt_linear[1] = tmp1; CCt_linear[2] = tmp1; CCt_linear[3] = tmp2;
	dgetrf(&nele1,&nele1,CCt_linear,&nele1,pivotArray1,&errinfo); 
	dgetri(&nele1,CCt_linear,&nele1,pivotArray1,lapackWorkspace1,&lwork1,&errinfo);
	#endif
	
	/*printf("\n");
	for(i=0;i<M;i++){
	for(j=0;j<M;j++){		
		printf("%g ", CCt[i*M+j]);
	}
	printf("\n");
	}*/

}

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
