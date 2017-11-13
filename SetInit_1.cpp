//static double vt[4] = {-0.3399810435848562648026658,0.3399810435848562648026658,-0.8611363115940525752239465,0.8611363115940525752239465};
//static double wt[4] = {0.6521451548625461426269361,0.6521451548625461426269361,0.3478548451374538573730639,0.3478548451374538573730639};
double wt[5]={0.5688888888888889, 0.4786286704993665, 0.4786286704993665,0.2369268850561891, 0.2369268850561891};
double vt[5]={0., -0.5384693101056831,0.5384693101056831,-0.9061798459386640,0.9061798459386640};
	
double sinc(double x)
{
  double result;
  
  if(x==0.0)
    {
      result = 1.0;
    }
  else
    {
      result = sin(x)/x;
    }
  
  return result;
}

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

void trapezoidalRule(int nPoints, double *weight)
{
	int i;
	
	weight[0] = 0.5;
	weight[nPoints-1] = weight[0];
	
	for(i=1;i<nPoints-1;i++) weight[i] = 1.0;
}

double f_TS(double v1, double v2, double v3) //two-stream instability initial
{
  double r2=v1*v1+v2*v2+v3*v3;
  return r2*exp(-r2/2)/(2*PI*sqrt(2*PI))/3.;
}

double f_2Gauss(double v1, double v2, double v3)
{	
  double retn, sig=M_PI/10;
  retn = 0.5*(exp(-((v1-2*sig)*(v1-2*sig)+v2*v2+v3*v3)/(2*sig*sig))+exp(-((v1+2*sig)*(v1+2*sig)+v2*v2+v3*v3)/(2*sig*sig)))/(2*M_PI*sig*sig*sqrt(2*M_PI*sig*sig));
  return retn;
}

double Mw(double v1, double v2, double v3)
{
  double r2, T, retn;
  T=0.9; // when testing the nonlinear damping, T was chosen too small that the "effective grid" is  not fine enough 
  r2=v1*v1+v2*v2+v3*v3;
  retn = exp(-r2/(2*T))/(2*PI*T*sqrt(2*T*PI));
 // if(v1<=dv && v1>=-dv)retn=0.5/(dv*Lv*Lv);
  //else retn=0.;
  return retn;
}

void SetInit(double *U)
{
    int i, j1, j2, j3, k, m1,m2,m3,nt=5;
    double a=A_amp, c=k_wave;
    double tp, tp0, tp5, tmp0, tmp1, tmp2, tmp3, tmp4;	
    //#pragma omp parallel for private(k,j1,j2,j3,i,tmp0, tmp1, tmp2, tmp3, tmp4, tp0, tp5, tp) shared(U)
    for(j1=0;j1<Nv;j1++){
      for(j2=0;j2<Nv;j2++){
	for(j3=0;j3<Nv;j3++){
	  tmp0=0.; tmp1=0.; tmp2=0.; tmp3=0.; tmp4=0.;
	  for(m1=0;m1<nt;m1++){
	    for(m2=0;m2<nt;m2++){
	      for(m3=0;m3<nt;m3++){
		#ifdef Damping
		tp = wt[m1]*wt[m2]*wt[m3]*Mw(Gridv((double)j1)+0.5*dv*vt[m1], Gridv((double)j2)+0.5*dv*vt[m2], Gridv((double)j3)+0.5*dv*vt[m3]);
		#endif
		
		#ifdef TwoStream 
		tp = wt[m1]*wt[m2]*wt[m3]*f_2Gauss(Gridv((double)j1)+0.5*dv*vt[m1], Gridv((double)j2)+0.5*dv*vt[m2], Gridv((double)j3)+0.5*dv*vt[m3]);
		#endif
		
		tmp0 += tp;
		tmp1 += tp*0.5*vt[m1];
		tmp2 += tp*0.5*vt[m2];
		tmp3 += tp*0.5*vt[m3];
		tmp4 += tp*0.25*(vt[m1]*vt[m1] + vt[m2]*vt[m2]+ vt[m3]*vt[m3]);
	      }
	    }
	  }
	  tmp0 = tmp0*0.5*0.5*0.5; tmp1 = tmp1*0.5*0.5*0.5; tmp2 = tmp2*0.5*0.5*0.5; tmp3 = tmp3*0.5*0.5*0.5; tmp4 = tmp4*0.5*0.5*0.5;
	  for(i=0;i<Nx;i++){
	      k=i*size_v + (j1*Nv*Nv + j2*Nv + j3);	
	      tp0 = (dx + (sin(c*Gridx(i+0.5)) - sin(c*Gridx(i-0.5)))*a/c)*tmp0/dx;
	      tp5 = (dx + (sin(c*Gridx(i+0.5)) - sin(c*Gridx(i-0.5)))*a/c)*tmp4/dx;
	      U[k*6+0] = 19*tp0/4. - 15*tp5;
	      U[k*6+5] = 60*tp5 - 15*tp0;

	      U[k*6+1] = (0.5*(sin(c*Gridx(i+0.5)) + sin(c*Gridx(i-0.5))) + (cos(c*Gridx(i+0.5)) - cos(c*Gridx(i-0.5)))/(c*dx))*(a/c)*tmp0*12./dx;
	      U[k*6+2] = (dx + (sin(c*Gridx(i+0.5)) - sin(c*Gridx(i-0.5)))*a/c)*tmp1*12/dx;
	      U[k*6+3] = (dx + (sin(c*Gridx(i+0.5)) - sin(c*Gridx(i-0.5)))*a/c)*tmp2*12/dx;
	      U[k*6+4] = (dx + (sin(c*Gridx(i+0.5)) - sin(c*Gridx(i-0.5)))*a/c)*tmp3*12/dx; 
	  }
	}
      }
    }
}

#ifdef MPI
void setInit_spectral(double *U, double **f)
{
  int i, j1, j2, j3, k, l, m ,n;
  for(i=chunk_Nx*myrank_mpi;i<chunk_Nx*(myrank_mpi+1) && i<Nx;i++){  
    for(l=0;l<N;l++){
      j1 = (l*h_v)/dv; // integer part = floor() for non-negative integers.
      if(j1==Nv)j1=Nv-1; // let the right end point lie in the last element
      for(m=0;m<N;m++){
		j2 = (m*h_v)/dv;
		if(j2==Nv)j2=Nv-1;
		for(n=0;n<N;n++){
			j3 = (n*h_v)/dv;
			if(j3==Nv)j3=Nv-1;
			k=i*size_v + (j1*Nv*Nv + j2*Nv + j3); // determine in which element the Fourier nodes lie	  
			f[i%chunk_Nx][l*N*N+m*N+n] = U[k*6+0] + U[k*6+2]*(v[l]-Gridv((double)j1))/dv + U[k*6+3]*(v[m]-Gridv((double)j2))/dv + U[k*6+4]*(v[n]-Gridv((double)j3))/dv + U[k*6+5]*( ((v[l]-Gridv((double)j1))/dv)*((v[l]-Gridv((double)j1))/dv) + ((v[m]-Gridv((double)j2))/dv)*((v[m]-Gridv((double)j2))/dv) + ((v[n]-Gridv((double)j3))/dv)*((v[n]-Gridv((double)j3))/dv) ); 
		  //BUG: index was "l*N*N+m*N+n*N" !!!!!!
		}
      }
    }
  }   
}
#else
void setInit_spectral(double *U, double **f)
{
  int i, j1, j2, j3, k, l, m ,n;  
  for(l=0;l<N;l++){
      j1 = (l*h_v)/dv; // integer part = floor() for non-negative integers.
      if(j1==Nv)j1=Nv-1; // let the right end point lie in the last element
      for(m=0;m<N;m++){
	  j2 = (m*h_v)/dv;
	  if(j2==Nv)j2=Nv-1;
		for(n=0;n<N;n++){
		  j3 = (n*h_v)/dv;
		  if(j3==Nv)j3=Nv-1;
		  for(i=0;i<Nx;i++){
		  k=i*size_v + (j1*Nv*Nv + j2*Nv + j3); // determine in which element the Fourier nodes lie	  
		  f[i][l*N*N+m*N+n] = U[k*6+0] + U[k*6+2]*(v[l]-Gridv((double)j1))/dv + U[k*6+3]*(v[m]-Gridv((double)j2))/dv + U[k*6+4]*(v[n]-Gridv((double)j3))/dv + U[k*6+5]*( ((v[l]-Gridv((double)j1))/dv)*((v[l]-Gridv((double)j1))/dv) + ((v[m]-Gridv((double)j2))/dv)*((v[m]-Gridv((double)j2))/dv) + ((v[n]-Gridv((double)j3))/dv)*((v[n]-Gridv((double)j3))/dv) ); 
		  //BUG: index was "l*N*N+m*N+n*N" !!!!!!
		  }
        }
      }
    }   
}
#endif
