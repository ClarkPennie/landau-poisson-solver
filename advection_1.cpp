//advection and moments routines

double wt[5]={0.5688888888888889, 0.4786286704993665, 0.4786286704993665,0.2369268850561891, 0.2369268850561891};				// weights for Gaussian quadrature
double vt[5]={0., -0.5384693101056831,0.5384693101056831,-0.9061798459386640,0.9061798459386640};								// node values for Gaussian quadrature over the interval [-1,1]

double Gridv(double m){ //v in [-Lv,Lv]
	return (-Lv+(m+0.5)*dv);
}

double Gridx(double m){ // x in [0,Lx]  (returns the x value at the mth discrete space-step, in the middle of the cell I_m)
	return (m+0.5)*dx;
}

double rho_x(double x, double *U, int i) // for x in I_i
{
  int j, k;
  double tmp=0.;
  //#pragma omp parallel for shared(U) reduction(+:tmp)
  for(j=0;j<size_v;j++){
	k=i*size_v + j;
	tmp += U[k*6+0] + U[k*6+1]*(x-Gridx((double)i))/dx + U[k*6+5]/4.;
  }
  
  return dv*dv*dv*tmp;
}

double rho(double *U, int i) //  \int f(x,v) dxdv in I_i * K_j
{
  int j, k;
  double tmp=0.;
 // #pragma omp parallel for shared(U) reduction(+:tmp)
  for(j=0;j<size_v;j++){
	k=i*size_v + j;
	tmp += U[k*6+0] + U[k*6+5]/4.;
  }
  
  return dx*dv*dv*dv*tmp;
}

double computePhi_x_0(double *U) // compute the constant coefficient of x in phi, which is actually phi_x(0) (Calculate C_E in the paper -between eq. 52 & 53?)
{
	int i, j, k, m, q;
	double tmp=0.;

	//#pragma omp parallel for private(j,q,m,k) shared(U) reduction(+:tmp) //reduction may change the final result a little bit
	for(j=0;j<size_v;j++){
		//j = j1*Nv*Nv + j2*Nv + j3;
		for(q=0;q<Nx;q++){
			for(m=0;m<q;m++){
				k=m*size_v + j;
				tmp += U[k*6] + U[k*6+5]/4.;
			}
			k=q*size_v + j;
			tmp += 0.5*(U[k*6+0] + U[k*6+5]/4.) - U[k*6+1]/12.;
		}
	}
	tmp = tmp*scalev*dx*dx;

	return 0.5*Lx - tmp/Lx;
}

double computePhi(double *U, double x, int ix)											// function to compute the potential Phi at a position x, contained in [x_(ix-1/2), x_(ix+1/2)]
{
	int i_out, i, j1, j2, j3, iNNN, j1NN, j2N, k;										// declare counters i_out (for the outer sum of i values), i, j1, j2, j3 for summing the contribution from cell I_i x K_(j1, j2, j3), iNNN (the value of i*Nv^3), j1NN (the value of j1*Nv^2), j2N (the value of j2*N) & k (the location in U of the cell I_i x K_(j1, j2, j3))
	double retn, sum1, sum3, sum4, x_diff, x_diff_mid, x_diff_sq, x_eval, C_E;			// declare retn (the value of Phi returned at the end), sum1 (the value of the first two sums), sum3 (the value of the third sum), sum4 (the value of the fourth sum), x_diff (the value of x - x_(ix-1/2)), x_diff_mid (the value of x - x_ix), x_diff_sq (the value of x_diff^2), x_eval (the value associated to the integral of (x - x_i)^2) & C_E (the value of the constant in the formula for phi)
	sum1 = 0;
	sum3 = 0;
	sum4 = 0;
	retn = 0;
	x_diff = x - Gridx(ix-0.5);
	x_diff_mid = x - Gridx(ix);
	x_diff_sq = x_diff*x_diff;
	x_eval = x_diff_mid*x_diff_mid*x_diff_mid/(6.*dx) - dx*x_diff_mid/8. - dx*dx/24.;

	for(i_out = 0; i_out < ix; i_out++)
	{
		for(i = 0; i < i_out; i++)
		{
			iNNN = i*Nv*Nv*Nv;
			for(j1 = 0; j1 < Nv; j1++)
			{
				j1NN = j1*Nv*Nv;
				for(j2 = 0; j2 < Nv; j2++)
				{
					j2N = j2*Nv;
					for(j3 = 0; j3 < Nv; j3++)
					{
						k = iNNN + j1NN + j2N + j3;
						sum1 += U[6*k] + U[6*k+5]/4.;
					}
				}
			}
		}
		iNNN = i_out*Nv*Nv*Nv;
		for(j1 = 0; j1 < Nv; j1++)
		{
			j1NN = j1*Nv*Nv;
			for(j2 = 0; j2 < Nv; j2++)
			{
				j2N = j2*Nv;
				for(j3 = 0; j3 < Nv; j3++)
				{
					k = iNNN + j1NN + j2N + j3;
					sum1 += U[6*k]/2 - U[6*k+1]/12. +  U[6*k+5]/8;
				}
			}
		}
	}
	sum1 = sum1*dx*dx;
	for(i = 0; i < i_out; i++)
	{
		iNNN = i*Nv*Nv*Nv;
		for(j1 = 0; j1 < Nv; j1++)
		{
			j1NN = j1*Nv*Nv;
			for(j2 = 0; j2 < Nv; j2++)
			{
				j2N = j2*Nv;
				for(j3 = 0; j3 < Nv; j3++)
				{
					k = iNNN + j1NN + j2N + j3;
					sum3 += (U[6*k] + U[6*k+5]/4.);
				}
			}
		}
	}
	sum3 = sum3*dx*x_diff;
	iNNN = ix*Nv*Nv*Nv;
	for(j1 = 0; j1 < Nv; j1++)
	{
		j1NN = j1*Nv*Nv;
		for(j2 = 0; j2 < Nv; j2++)
		{
			j2N = j2*Nv;
			for(j3 = 0; j3 < Nv; j3++)
			{
				k = iNNN + j1NN + j2N + j3;
				sum4 += U[6*k]*x_diff_sq/2. + U[6*k+1]*x_eval +  U[6*k+5]*x_diff_sq/8.;
			}
		}
	}

	C_E = computePhi_x_0(U);
	retn = (sum1 + sum3 + sum4)*dv*dv*dv - x*x/2 - C_E*x;
	return retn;
}

void PrintPhiVals(double *U, FILE *phifile)														// function to print the values of the potential and the density in the file tagged as phifile at the given timestep
{
	int i, np, nx, nv;																			// declare i (the index of the space cell), np (the number of points to evaluate in a given space/velocity cell), nx (a counter for the points in the space cell) & nv (a counter for the points in the velocity cell)
	double x_0, x_val, phi_val, rho_val, M_0, ddx;												// declare x_0 (the x value at the left edge of a given cell), x_val (the x value to be evaluated at), phi_val (the value of phi evaluated at x_val), rho_val (the value of rho evaluated at x_val) & ddx (the space between x values)

	np = 4;																						// set np to 4
	ddx = dx/np;																				// set ddx to the space cell width divided by np
	for(i=0; i<Nx; i++)
	{
		x_0 = Gridx((double)i - 0.5);															// set x_0 to the value of x at the left edge of the i-th space cell
		for (nx=0; nx<np; nx++)
		{
			x_val = x_0 + nx*ddx;																// set x_val to x_0 plus nx increments of width ddx

			phi_val = computePhi(U, x_val, i);													// calculate the value of phi, evaluated at x_val by using the function in the space cell
			rho_val = rho_x(x_val, U, i);														// calculate the value of rho, evaluated at x_val by using the function in the space cell
			M_0 = rho_val/(sqrt(1.8*PI));
			fprintf(phifile, "%11.8g %11.8g %11.8g \n", phi_val, rho_val, M_0);								// in the file tagged as phifile, print the value of the potential phi(t, x) and the density rho(t, x)
		}
	}
}

double computeC_rho(double *U, int i) // sum_m=0..i-1 int_{I_m} rho(z)dz   (Calculate integral of rho_h(z,t) from 0 to x as in eq. 53)
{
	double retn=0.;
	int j, k, m;
	for (m=0;m<i;m++){ //BUG: was "m < i-1"
		for(j=0;j<size_v;j++){
			k = m*size_v + j;
			retn += U[k*6+0] + U[k*6+5]/4.;
		}
	}
	
	retn *= dx*scalev;
	return retn;
}

double Int_Int_rho(double *U, int i) // \int_{I_i} [ \int^{x}_{x_i-0.5} rho(z)dz ] dx
{
  int j, k;
  double retn=0.;
  for(j=0;j<size_v;j++){
    k=i*size_v + j;
    retn += 0.5*(U[k*6+0] + U[k*6+5]/4.) - U[k*6+1]/12.;
  }
  
  return retn*dx*dx*scalev;  
}

double Int_Int_rho1st(double *U, int i)// \int_{I_i} [(x-x_i)/delta_x * \int^{x}_{x_i-0.5} rho(z)dz ] dx
{
 int j, k;
 double retn=0.;
 for(j=0;j<size_v;j++){
    k=i*size_v + j;
    retn += (U[k*6+0] + U[k*6+5]/4.)/12.;
  }
  return retn*dx*dx*scalev; 
}

/*double Int_Cumulativerho(double **U, int i)// \int_{I_i} [ \int^{x}_{0} rho(z)dz ] dx
{
  double retn=0., cp, tmp;
  cp = computeC_rho(U,i);
  tmp = Int_Int_rho(U,i);
  
  retn = dx*cp + tmp;
  return retn;  
}

double Int_Cumulativerho_sqr(double **U, int i)// \int_{I_i} [ \int^{x}_{0} rho(z)dz ]^2 dx
{
  double retn=0., cp, tmp1, tmp2, tmp3, c1=0., c2=0.;
  int j, k;
  cp = computeC_rho(U,i);
  tmp1 = cp*cp*dx;
  tmp2 = 2*cp*Int_Int_rho(U,i);
  for(j=0;j<size_v;j++){
    k=i*size_v + j;
    c1 += U[k][0] + U[k][5]/4.;
    c2 += U[k][1];
  }
  c2 *= dx/2.;
  tmp3 = pow(dv, 6)* ( c1*c1*dx*dx*dx/3. + c2*c2*dx/30. + c1*c2*(dx*Gridx((double)i)/6. - dx*dx/4.) );
  retn = tmp1 + tmp2 + tmp3;
  return retn;  
}*/

double Int_E(double *U, int i) // \int_i E dx      // Function to calculate the integral of E_h w.r.t. x over the interval I_i = [x_(i-1/2), x_(i+1/2))
{
	int m, j, k;
	double tmp=0., result;
	//#pragma omp parallel for shared(U) reduction(+:tmp)
	for(j=0;j<size_v;j++){
		for(m=0;m<i;m++){	
			k=m*size_v + j;
			tmp += U[k*6+0] + U[k*6+5]/4.;
		}
		k=i*size_v + j;
		tmp += 0.5*(U[k*6+0] + U[k*6+5]/4.) - U[k*6+1]/12.;		
	}

	//ce = computePhi_x_0(U);
	result = -ce*dx - tmp*dx*dx*scalev + Gridx((double)i)*dx;

	return result;
}

double Int_E1st(double *U, int i) // \int_i E*(x-x_i)/delta_x dx
{
	int j, k;
	double tmp=0., result;
	//#pragma omp parallel for reduction(+:tmp)
	for(j=0;j<size_v;j++){
		k=i*size_v + j;
		tmp += U[k*6+0] + U[k*6+5]/4.;
	}
	tmp = tmp*scalev;
	
	result = (1-tmp)*dx*dx/12.;

	return result;
}

double Int_fE(double *U, int i, int j) // \int f * E(f) dxdv on element I_i * K_j
{
	double retn=0.;
	int k;
	k = i*size_v + j;

	//retn = (U[k][0] + U[k][5]/4.)*Int_E(U,i) + U[k][1]*Int_E1st(U,i);
	retn = (U[k*6+0] + U[k*6+5]/4.)*intE[i] + U[k*6+1]*intE1[i];
	return retn*scalev;
}



double Int_E2nd(double *U, int i) // \int_i E* [(x-x_i)/delta_x]^2 dx
{
    int m, j, j1, j2, j3, k;
    double c1=0., c2=0., result;
  
    //cp = computeC_rho(U,i); ce = computePhi_x_0(U);

    for(j=0;j<size_v;j++){
	    k=i*size_v + j;
	    c1 += U[k*6+0] + U[k*6+5]/4.;
	    c2 += U[k*6+1];
    }
    c2 *= dx/2.;				
    
    result = (-cp[i] - ce+scalev*(c1*Gridx(i-0.5) + 0.25*c2))*dx/12. + (1-scalev*c1)*dx*Gridx((double)i)/12. - scalev*c2*dx/80.; //BUG: missed -cp

    return result;
}

double computeMass(double *U)
{
  int k;
  double tmp=0.;
  #pragma omp parallel for shared(U) reduction(+:tmp)
  for(k=0;k<Nx*size_v;k++) tmp += U[k*6+0] + U[k*6+5]/4.;
  
  return tmp*dx*scalev;
}

void computeMomentum(double *U, double *a)
{
  int k, i, j,j1,j2,j3; //k=i*Nv*Nv*Nv + (j1*Nv*Nv + j2*Nv + j3);	
  double tmp1=0., tmp2=0., tmp3=0.; 
  a[0]=0.; a[1]=0.; a[2]=0.; // the three momentum
  #pragma omp parallel for shared(U) reduction(+:tmp1, tmp2, tmp3)  //reduction directive may change the result a little bit
  for(k=0;k<Nx*size_v;k++){
    j=k%size_v; i=(k-j)/size_v;
    j3=j%Nv; j2=((j-j3)%(Nv*Nv))/Nv; j1=(j-j3-j2*Nv)/(Nv*Nv);
    tmp1 += Gridv((double)j1)*dv*U[k*6+0] + U[k*6+2]*dv*dv/12. + U[k*6+5]*Gridv((double)j1)*dv/4.;
    tmp2 += Gridv((double)j2)*dv*U[k*6+0] + U[k*6+3]*dv*dv/12. + U[k*6+5]*Gridv((double)j2)*dv/4.;
    tmp3 += Gridv((double)j3)*dv*U[k*6+0] + U[k*6+4]*dv*dv/12. + U[k*6+5]*Gridv((double)j3)*dv/4.;
  }
  a[0]=tmp1*dx*dv*dv; a[1]=tmp2*dx*dv*dv; a[2]=tmp3*dx*dv*dv; 
}

double computeKiE(double *U)
{
  int k, i, j,j1,j2,j3; 
  double tmp=0., tp=0., tp1=0.;  
  //#pragma omp parallel for private(k,i,j,j1,j2,j3,tp, tp1) shared(U) reduction(+:tmp)
  for(k=0;k<Nx*size_v;k++){
    j=k%size_v; i=(k-j)/size_v;
    j3=j%Nv; j2=((j-j3)%(Nv*Nv))/Nv; j1=(j-j3-j2*Nv)/(Nv*Nv);
    //tp = ( pow(Gridv(j1+0.5), 3)- pow(Gridv(j1-0.5), 3) + pow(Gridv(j2+0.5), 3)- pow(Gridv(j2-0.5), 3) + pow(Gridv(j3+0.5), 3)- pow(Gridv(j3-0.5), 3) )/3.;
    tp1 = Gridv((double)j1)*Gridv((double)j1) + Gridv((double)j2)*Gridv((double)j2) + Gridv((double)j3)*Gridv((double)j3);
    //tmp += U[k][0]*tp + (Gridv(j1)*U[k][2]+Gridv(j2)*U[k][3]+Gridv(j3)*U[k][4])*dv*dv/6. + U[k][5]*( dv*(dv*dv*3./80. + tp1/12.) + tp/6. );	 
    tmp += U[k*6+0]*(tp1 + dv*dv/4.)*dv + (Gridv(j1)*U[k*6+2]+Gridv(j2)*U[k*6+3]+Gridv(j3)*U[k*6+4])*dv*dv/6. + U[k*6+5]*( dv*dv*dv*19./240. + tp1*dv/4.);	 
  }
  tmp *= dx*dv*dv; 
  return 0.5*tmp;
}

double computeKiEratio(double *U, int *NegVals)
{
	int k, i, j,j1,j2,j3;
	double KiEpos, KiEneg;																	// declare KiEpos (the kinetic energy where f is positive all over the cells) & KiEneg (the kinetic energy where f goes negative in the cells)
	double tmp=0., tp=0., tp1=0.;
	KiEpos = 0;
	KiEneg = 0;
	#pragma omp parallel for private(k,i,j,j1,j2,j3,tp, tp1) shared(U) reduction(+:KiEpos,KiEneg)
	for(k=0;k<Nx*size_v;k++)
	{
		if(NegVals[k] == 0)
		{
			j=k%size_v; i=(k-j)/size_v;
			j3=j%Nv; j2=((j-j3)%(Nv*Nv))/Nv; j1=(j-j3-j2*Nv)/(Nv*Nv);
			tp1 = Gridv((double)j1)*Gridv((double)j1) + Gridv((double)j2)*Gridv((double)j2) + Gridv((double)j3)*Gridv((double)j3);
			KiEpos += U[k*6+0]*(tp1 + dv*dv/4.)*dv + (Gridv(j1)*U[k*6+2]+Gridv(j2)*U[k*6+3]+Gridv(j3)*U[k*6+4])*dv*dv/6. + U[k*6+5]*( dv*dv*dv*19./240. + tp1*dv/4.);
		}
		else
		{
			j=k%size_v; i=(k-j)/size_v;
			j3=j%Nv; j2=((j-j3)%(Nv*Nv))/Nv; j1=(j-j3-j2*Nv)/(Nv*Nv);
			tp1 = Gridv((double)j1)*Gridv((double)j1) + Gridv((double)j2)*Gridv((double)j2) + Gridv((double)j3)*Gridv((double)j3);
			KiEneg += U[k*6+0]*(tp1 + dv*dv/4.)*dv + (Gridv(j1)*U[k*6+2]+Gridv(j2)*U[k*6+3]+Gridv(j3)*U[k*6+4])*dv*dv/6. + U[k*6+5]*( dv*dv*dv*19./240. + tp1*dv/4.);
		}
	}
	return KiEneg/KiEpos;
}

double computeEleE(double *U)
{
  int k, i, j;
  double retn, tmp1=0., tmp2=0., tmp3=0., tmp4=0., tmp5=0., tmp6=0., tmp7=0., tp1, tp2, c;
  double ce1, cp1;
  ce1 = computePhi_x_0(U);
  
  tmp1 = ce1*ce1*Lx;
  tmp2 = Lx*Lx*Lx/3.; tmp3 = -ce1*Lx*Lx;
  
  //#pragma omp parallel for private(j,k, i, tp1, tp2, cp1, c) shared(U) reduction(+:tmp4, tmp5, tmp6)
  for(i=0;i<Nx;i++){
    c = Int_Int_rho(U,i);
    cp1 = computeC_rho(U,i);
    tmp4 += dx*cp1 + c;   
    tmp5 += dx*Gridx((double)i)*cp1;
    tp1=0.; tp2=0.;
    for(j=0;j<size_v;j++){
      k = i*size_v + j;
      tp1 += (U[k*6+0] + U[k*6+5]/4.);
      tp2 += U[k*6+1];
    }
    tmp5 += scalev* (tp1*( (pow(Gridx(i+0.5), 3) - pow(Gridx(i-0.5), 3))/3. - Gridx(i-0.5)*Gridx((double)i)*dx ) - tp2 * dx*dx*Gridx((double)i)/12.);
    
    tp2 *= dx/2.;
    tmp6 +=  cp1*cp1*dx + 2*cp1*c + pow(dv, 6)* ( tp1*tp1*dx*dx*dx/3. + tp2*tp2*dx/30. - tp1*tp2*dx*dx/6.);//+ tp1*tp2*(dx*Gridx((double)i)/6. - dx*dx/4.) ); //Int_Cumulativerho_sqr(i);
  }
  retn = tmp1 + tmp2 + tmp3 + 2*ce1*tmp4 - 2*tmp5 + tmp6;
  return 0.5*retn;
}

double computeEntropy(double *U)																									// function to compute the entropy for the given approximate solution, namely \int f log(f) dx dv, over Omega_x x Omega_v
{
	int i, j1, j2, j3, iNNN, j1NN, j2N, k, nx, nv1, nv2, nv3;																		// declare i (the index of the space cell), j1, j2, j3 (the indices of the velocity cell), iNNN (to store i*Nv^3), j1NN (to store j1*Nv^2), j2N (to store j2*Nv), k (the location of the given cell in U), nx (a counter for the space cell) & nv1, nv2, nv3 (counters for the velocity cell)
	double x_0, v1_0, v2_0, v3_0, x_val, v1_val, v2_val, v3_val, f_val;																// declare x_0, v1_0, v2_0, v3_0 (to store the coordinates in the middle of the current cell), x_val, v1_val, v2_val, v3_val (to store the x & v values to be evaluated at) & f_val (to store the value of the function evaluated at the current x & v values)
	double Ent;																														// declare Ent (the value of the entropy)
	Ent = 0;																														// set Ent to 0 before any quadrature has began
	#pragma omp parallel for private(k,i,j1,j2,j3,iNNN,j1NN,j2N,nx,nv1,nv2,nv3,x_0,v1_0,v2_0,v3_0,x_val,v1_val,v2_val,v3_val,f_val) shared(U,size_v,Nv,dv,dx,vt,wt) reduction(+:Ent)
	for(i=0;i<Nx;i++)																												// loop through the space cells
	{
		iNNN = i*size_v;																											// set iNNN to i*Nv^3
		x_0 = Gridx((double)i);																												// set x_0 to the value of x at the center of the ith space cell

		for(j1=0;j1<Nv;j1++)																										// loop through the velocity cells in the v1 direction
		{
			j1NN = j1*Nv*Nv;																										// set j1NN to j1*Nv^2
			v1_0 = Gridv((double)j1);																										// set v1_0 to the value of v1 at the center of the j1th velocity cell in that direction
			for(j2=0;j2<Nv;j2++)																									// loop through the velocity cells in the v2 direction
			{
				j2N = j2*Nv;																										// set j2N to j2*Nv
				v2_0 = Gridv((double)j2);																									// set v2_0 to the value of v2 at the center of the j2th velocity cell in that direction
				for(j3=0;j3<Nv;j3++)																								// loop through the velocity cells in the v3 direction
				{
					k = iNNN + j1NN + j2N + j3;																						// set k to i*Nv^3 + j1*Nv^2 + j2*Nv + j3
					v3_0 = Gridv((double)j3);																								// set v3_0 to the value of v3 at the center of the j3th velocity cell in that direction
					for(nx=0;nx<5;nx++)																								// loop through the five quadrature points in the x direction of the cell
					{
						x_val = x_0 + 0.5*vt[nx]*dx;																					// set x_val to the nx-th quadrature point in the cell
						for(nv1=0;nv1<5;nv1++)																						// loop through the five quadrature points in the v1 direction of the cell
						{
							v1_val = v1_0 + 0.5*vt[nv1]*dv;																			// set x_val to the nv1-th quadrature point in the cell
							for(nv2=0;nv2<5;nv2++)																					// loop through the five quadrature points in the v2 direction of the cell
							{
								v2_val = v2_0 + 0.5*vt[nv2]*dv;																		// set x_val to the nv2-th quadrature point in the cell
								for(nv3=0;nv3<5;nv3++)																				// loop through the five quadrature points in the v3 direction of the cell
								{
									v3_val = v3_0 + 0.5*vt[nv3]*dv;																	// set x_val to the nv3-th quadrature point in the cell
									f_val = U[k*6+0] + U[k*6+1]*(x_val-x_0)/dx + U[k*6+2]*(v1_val-v1_0)/dv +
											U[k*6+3]*(v2_val-v2_0)/dv + U[k*6+4]*(v3_val-v3_0)/dv +
											U[k*6+5]*(((v1_val-v1_0)/dv)*((v1_val-v1_0)/dv)+((v2_val-v2_0)/dv)*((v2_val-v2_0)/dv)
													+((v3_val-v3_0)/dv)*((v3_val-v3_0)/dv));										// set f_val to the evaluation of the approximation at (x_val,v1_val,v2_val,v3_val)
									if(f_val > 0)																					// only do this if f(x_val,v1_val,v2_val,v3_val) > 0 so that the log can be evaluated
									{
										Ent += wt[nx]*wt[nv1]*wt[nv2]*wt[nv3]*f_val*log(f_val);										// add f(x_val,v1_val,v2_val,v3_val)*log(f(x_val,v1_val,v2_val,v3_val)) times the quadrature weights corresponding to the current x & v values to the current quadrature value
									}
								}
							}
						}
					}
				}
			}
		}
	}
	Ent = Ent*0.5*dv*0.5*dv*0.5*dv*0.5*dx;																							// multiply the quadrature result by 0.5*dv*0.5*dv*0.5*dv*0.5*dx, since each quadrature is done over intervals of width dv (three times) or dx (once) instead of the standard interval [-1,1]
	return Ent;
}

double computeEntropy_wAvg(double *AvgVals)																								// function to compute the entropy for the given approximate solution, namely \int f log(f) dx dv, over Omega_x x Omega_v, with the average values of f on the cells
{
	int i, j1, j2, j3, iNNN, j1NN, j2N, k;																							// declare i (the index of the space cell), j1, j2, j3 (the indices of the velocity cell), iNNN (to store i*Nv^3), j1NN (to store j1*Nv^2), j2N (to store j2*Nv) & (the location of the given cell in U)
	//double x_0, v1_0, v2_0, v3_0, x_val, v1_val, v2_val, v3_val, f_val;																// declare x_0, v1_0, v2_0, v3_0 (to store the coordinates in the middle of the current cell), x_val, v1_val, v2_val, v3_val (to store the x & v values to be evaluated at) & f_val (to store the value of the function evaluated at the current x & v values)
	double f_val, Ent;																														// declare Ent (the value of the entropy)
	Ent = 0;																														// set Ent to 0 before any quadrature has began
	#pragma omp parallel for private(k,i,j1,j2,j3,iNNN,j1NN,j2N,f_val) shared(AvgVals,size_v,Nx,Nv) reduction(+:Ent)
	for(i=0;i<Nx;i++)																												// loop through the space cells
	{
		iNNN = i*size_v;																											// set iNNN to i*Nv^3
		for(j1=0;j1<Nv;j1++)																										// loop through the velocity cells in the v1 direction
		{
			j1NN = j1*Nv*Nv;																										// set j1NN to j1*Nv^2
			for(j2=0;j2<Nv;j2++)																									// loop through the velocity cells in the v2 direction
			{
				j2N = j2*Nv;																										// set j2N to j2*Nv
				for(j3=0;j3<Nv;j3++)																								// loop through the velocity cells in the v3 direction
				{
					k = iNNN + j1NN + j2N + j3;																						// set k to i*Nv^3 + j1*Nv^2 + j2*Nv + j3
					f_val = fabs(AvgVals[k]);																						// set f_val to the absolute value of the average value of f on the cell, stored in AvgVals[k]
					Ent += f_val*log(f_val);																						// add f(x_val,v1_val,v2_val,v3_val)*log(f(x_val,v1_val,v2_val,v3_val)
				}
			}
		}
	}
	Ent = Ent*dv*dv*dv*dx;																											// multiply the result by dv*dv*dv*dx, since each addition to the integral should be multiplied by the volume of the cell
	return Ent;
}

double computeRelEntropy(double *U, double *f_equi)																									// function to compute the entropy for the given approximate solution, namely \int f log(f) dx dv, over Omega_x x Omega_v
{
	int i, j1, j2, j3, iNNN, j1NN, j2N, k, k_equi, nx, nv1, nv2, nv3;																// declare i (the index of the space cell), j1, j2, j3 (the indices of the velocity cell), iNNN (to store i*Nv^3), j1NN (to store j1*Nv^2), j2N (to store j2*Nv), k (the location of the given cell in U), k_equi (the location of the current quadrature point in f_equi) &  nx, nv1, nv2, nv3 (counters for the Gaussian quadrature)
	double x_0, v1_0, v2_0, v3_0, x_val, v1_val, v2_val, v3_val, f_val;																// declare x_0, v1_0, v2_0, v3_0 (to store the coordinates in the middle of the current cell), x_val, v1_val, v2_val, v3_val (to store the x & v values to be evaluated at) & f_val (to store the value of the function evaluated at the current v values)
	double Ent;																														// declare Ent (the value of the entropy)
	Ent = 0;																														// set Ent to 0 before any quadrature has began
	#pragma omp parallel for private(i,j1,j2,j3,iNNN,j1NN,j2N,k,k_equi,nx,nv1,nv2,nv3,x_0,v1_0,v2_0,v3_0,x_val,v1_val,v2_val,v3_val,f_val) shared(U,f_equi,size_v,Nv,dv,vt,wt) reduction(+:Ent)
	for(i=0;i<Nx;i++)																												// loop through the space cells
	{
		iNNN = i*size_v;																		// set iNNN to i*Nv^3
		x_0 = Gridx((double)i);																	// set x_0 to the value of x at the center of the ith space cell
		for(j1=0;j1<Nv;j1++)																										// loop through the velocity cells in the v1 direction
		{
			j1NN = j1*Nv*Nv;																										// set j1NN to j1*Nv^2
			v1_0 = Gridv((double)j1);																										// set v1_0 to the value of v1 at the center of the j1th velocity cell in that direction
			for(j2=0;j2<Nv;j2++)																									// loop through the velocity cells in the v2 direction
			{
				j2N = j2*Nv;																										// set j2N to j2*Nv
				v2_0 = Gridv((double)j2);																									// set v2_0 to the value of v2 at the center of the j2th velocity cell in that direction
				for(j3=0;j3<Nv;j3++)																								// loop through the velocity cells in the v3 direction
				{
					k = iNNN + j1NN + j2N + j3;																							// set k to j1*Nv^2 + j2*Nv + j3
					v3_0 = Gridv((double)j3);																						// set v3_0 to the value of v3 at the center of the j3th velocity cell in that direction
					for(nx=0;nx<5;nx++)															// loop through the five quadrature points in the x direction of the cell
					{
						x_val = x_0 + 0.5*vt[nx]*dx;											// set x_val to the nx-th quadrature point in the cell
						for(nv1=0;nv1<5;nv1++)																						// loop through the five quadrature points in the v1 direction of the cell
						{
							v1_val = v1_0 + 0.5*vt[nv1]*dv;																			// set v1_val to the nv1-th quadrature point in the cell
							for(nv2=0;nv2<5;nv2++)																					// loop through the five quadrature points in the v2 direction of the cell
							{
								v2_val = v2_0 + 0.5*vt[nv2]*dv;																		// set v2_val to the nv2-th quadrature point in the cell
								for(nv3=0;nv3<5;nv3++)																				// loop through the five quadrature points in the v3 direction of the cell
								{
									v3_val = v3_0 + 0.5*vt[nv3]*dv;																	// set v3_val to the nv3-th quadrature point in the cell
									f_val = U[k*6+0] + U[k*6+1]*(x_val-x_0)/dx + U[k*6+2]*(v1_val-v1_0)/dv +
											U[k*6+3]*(v2_val-v2_0)/dv + U[k*6+4]*(v3_val-v3_0)/dv +
											U[k*6+5]*(((v1_val-v1_0)/dv)*((v1_val-v1_0)/dv)+((v2_val-v2_0)/dv)*((v2_val-v2_0)/dv)
													+((v3_val-v3_0)/dv)*((v3_val-v3_0)/dv));										// set f_val to the evaluation of the approximation at (x_val,v1_val,v2_val,v3_val)
									if(f_val > 0)																					// only do this if f(v1_val,v2_val,v3_val) > 0 so that the log can be evaluated
									{
										k_equi = 5*(5*(5*(5*k + nx) + nv1) + nv2) + nv3;			// set k_equi to the location that this value of f should be stored in f
										/*
										if(nx==2 && nv1==2 && nv2==2 && nv3 == 2)
										{
											if(isnan(log(f_val/f_equi[k_equi])) != 0)
											{
												printf("(%d,%d,%d,%d): %g %g %g %g\n", i, j1, j2, j3, f_val, f_equi[k_equi], f_val/f_equi[k_equi], log(f_val/f_equi[k_equi]));
											}
										}
										*/
										if(isnan(log(f_val/f_equi[k_equi])) == 0)
										{
											Ent += wt[nx]*wt[nv1]*wt[nv2]*wt[nv3]*f_val*log(f_val/f_equi[k_equi]);									// add f(v1_val,v2_val,v3_val)*log(f(v1_val,v2_val,v3_val)/f_eq(v1_val,v2_val,v3_val)) times the quadrature weights corresponding to the current x & v values to the current quadrature value
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	Ent = Ent*0.5*dv*0.5*dv*0.5*dv*0.5*dx;																							// multiply the quadrature result by 0.5*dv*0.5*dv*0.5*dv, since each quadrature is done over intervals of width dv (three times) instead of the standard interval [-1,1]
	return Ent;
}

double I1(double *U, int k, int l) // Calculate the first inegtral in H_(i,j), namely \int v1*f*phi_x dxdv
{
  double result;
  int i, j1, j2, j3; // k=i*Nv^3 + (j1*Nv*Nv + j2*Nv + j3)
  int j_mod = k%size_v;
  j3 = j_mod%Nv;
  j2 = ((j_mod-j3)%(Nv*Nv))/Nv;
  j1 = (j_mod-j3-j2*Nv)/(Nv*Nv);
  i = (k-j_mod)/size_v;

  if(l==1) result = dv*dv*dv*( Gridv((double)j1)*U[k*6+0] + dv*U[k*6+2]/12. + U[k*6+5]*Gridv((double)j1)/4.);
  else result=0.;
  
  return result;
}

double I2(double *U, int k, int l) // Calculate the fourth integral in H_(i,j), namely \int E*f*phi_v1 dxdv
{
  double result;
  int i, j1, j2, j3; // k=i*Nv^3 + (j1*Nv*Nv + j2*Nv + j3)
  int j = k%size_v;
  //j3 = j_mod%Nv;
  //j2 = ((j_mod-j3)%(Nv*Nv))/Nv;
  //j1 = (j_mod-j3-j2*Nv)/(Nv*Nv);
  i = (k-j)/size_v;

  if(l==2) result = Int_fE(U,i,j)/dv;
  else if(l==5) result = U[k*6+2]*dv*dv*intE[i]/6.; 
  else result = 0.;
  
  return result;
}

double I3(double *U, int k, int l) 																			// Calculate the difference of the second and third integrals in H_(i,j), namely \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2
{
	double result, ur, ul;																					// declare result (the result of the integral to be returned), ur (used in the evaluation of gh^+/- on the right space cell edge) & ul (used in the evaluation of gh^+/- on the left space cell edge)
	int i, j1, j2, j3, iil, iir, kkl, kkr; 																	// declare i (the space cell coordinate), j1, j2, j3 (the coordinates of the velocity cell), iil (the cell from which the flux is flowing on the left of cell i in space), iir (the cell from which the flux is flowing on the right of cell i in space), kkl (the global index of the cell with coordinate (iil, j1, j2, j3)) & kkr (the global index of the cell with coordinate (iir, j1, j2, j3))
	int j_mod = k%size_v;																					// declare and calculate j_mod (the remainder when k is divided by size_v = Nv^3 - used to help determine the values of i, j1, j2 & j3 from the value of k)
	j3 = j_mod%Nv;																							// calculate j3 for the given k
	j2 = ((j_mod-j3)%(Nv*Nv))/Nv;																			// calculate j2 for the given k
	j1 = (j_mod-j3-j2*Nv)/(Nv*Nv);																			// calculate j1 for the given k
	i = (k-j_mod)/size_v;																					// calculate i for the given k

	if(j1<Nv/2)																								// do this if j1 < Nv/2 (so that the velocity in the v1 direction is negative)
	{
		iir=i+1; iil=i; 																					// set iir to the value of i+1 and iil to the value of i (as here the flow of information is from right to left so that gh^+ must be used at the cell edges)
		if(iir==Nx)iir=0; //periodic bc																		// if iir = Nx (the maximum value that can be obtained, since i = 0,1,...,Nx-1) and so this cell is at the right boundary, requiring information from the non-existent cell with space index Nx, since there are periodic boundary conditions, set iir = 0 and use the cell with space index 0 (i.e. the cell at the left boundary)
		kkr=iir*size_v + j_mod; 																			// calculate the value of kkr for this value of iir
		kkl=k;																								// set kkl to k (since iil = i)
		ur = -U[kkr*6+1]; 																					// set ur to the negative of the coefficient of the basis function with shape l which is non-zero in the cell with global index kkr (which corresponds to the evaluation of gh^+ at the right boundary and -ve since v_1 < 0 in here)
		ul = -U[kkl*6+1];																					// set ul to the negative of the coefficient of the basis function with shape l which is non-zero in the cell with global index kkl	(which corresponds to the evaluation of gh^+ at the left boundary and -ve since phi < 0 here)
	}
	else																									// do this if j1 >= Nv/2 (so that the velocity in the v1 direction is non-negative)
	{
		iir=i; iil=i-1;																						// set iir to the value of i and iil to the value of i-1 (as here the flow of information is from left to right so that gh^- must be used at the cell edges)
		if(iil==-1)iil=Nx-1; // periodic bc																	// if iil = -1 (the minimum value that can be obtained, since i = 0,1,...,Nx-1) and so this cell is at the left boundary, requiring information from the non-existent cell with space index -1, since there are periodic boundary conditions, set iil = Nx-1 and use the cell with space index Nx-1 (i.e. the cell at the right boundary)
		kkr=k; 																								// set kkr to k (since iir = i)
		kkl=iil*size_v + j_mod; 																			// calculate the value of kkl for this value of iil
		ur = U[kkr*6+1];																					// set ur to the value of the coefficient of the basis function with shape l which is non-zero in the cell with global index kkr (which corresponds to the evaluation of gh^- at the right boundary and +ve since v_1 >= 0 in here)
		ul = U[kkl*6+1];																					// set ul to the value of the coefficient of the basis function with shape l which is non-zero in the cell with global index kkl (which corresponds to the evalutaion of gh^- at the left boundary and +ve since v_r >= 0 in here)
	}
  
	if(l==0)result = dv*dv*dv*( (U[kkr*6+0]+0.5*ur - U[kkl*6+0]-0.5*ul)*Gridv((double)j1) + (U[kkr*6+2]-U[kkl*6+2])*dv/12. + (U[kkr*6+5]-U[kkl*6+5])*Gridv((double)j1)/4.);					// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 0 (i.e. constant) which is non-zero in the cell with global index k
	if(l==1)result = 0.5*dv*dv*dv*( (U[kkr*6+0]+0.5*ur + U[kkl*6+0]+0.5*ul)*Gridv((double)j1) + (U[kkr*6+2]+U[kkl*6+2])*dv/12. + (U[kkr*6+5]+U[kkl*6+5])*Gridv((double)j1)/4.);				// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 1 (i.e. linear in x) which is non-zero in the cell with global index k
	if(l==2)result = dv*dv*(( (U[kkr*6+0]-U[kkl*6+0])*dv*dv + (ur-ul)*0.5*dv*dv + (U[kkr*6+2]-U[kkl*6+2])*dv*Gridv((double)j1))/12. + (U[kkr*6+5]-U[kkl*6+5])*dv*dv*19./720.);				// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 2 (i.e. linear in v_1) which is non-zero in the cell with global index k
	if(l==3)result = (U[kkr*6+3]-U[kkl*6+3])*Gridv((double)j1)*dv*dv*dv/12.;																												// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 3 (i.e. linear in v_2) which is non-zero in the cell with global index k
	if(l==4)result = (U[kkr*6+4]-U[kkl*6+4])*Gridv((double)j1)*dv*dv*dv/12.;																												// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 4 (i.e. linear in v_3) which is non-zero in the cell with global index k
	if(l==5)result = dv*dv*dv*((U[kkr*6+0] + 0.5*ur - U[kkl*6+0]-0.5*ul)*Gridv((double)j1)/4. + (U[kkr*6+2]-U[kkl*6+2])*dv*19./720. + (U[kkr*6+5]-U[kkl*6+5])*Gridv((double)j1)*19./240.);	// calculate \int_j v1*gh*phi dv at interface x=x_i+1/2 - \int_j v1*gh*phi dv at interface x=x_i-1/2 for the basis function with shape 5 (i.e. modulus of v) which is non-zero in the cell with global index k

	return result;
}

double I5(double *U, int k, int l) 	// Calculate the difference of the fifth and sixth integrals in H_(i,j), namely \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2
{
	double result, ur, ul;																			// declare result (the result of the integral to be returned), ur (used in the evaluation of gh^+/- on the right velocity cell edge) & ul (used in the evaluation of gh^+/- on the left velocity cell edge)
	int i, j1, j2, j3, j1r, j1l, kkr, kkl; 															// declare i (the space cell coordinate), j1, j2, j3 (the coordinates of the velocity cell), j1l (the cell from which the flux is flowing on the left of the cell with coordinate j1 in the v1 direction), j1r (the cell from which the flux is flowing on the right of the cell with coordinate j1 in the v1 direction), kkl (the global index of the cell with coordinate (i, j1l, j2, j3)) & kkr (the global index of the cell with coordinate (i, j1r, j2, j3))
	int j_mod = k%size_v;																			// declare and calculate j_mod (the remainder when k is divided by size_v = Nv^3 - used to help determine the values of i, j1, j2 & j3 from the value of k)
	j3 = j_mod%Nv;																					// calculate j3 for the given k
	j2 = ((j_mod-j3)%(Nv*Nv))/Nv;																	// calculate j2 for the given k
	j1 = (j_mod-j3-j2*Nv)/(Nv*Nv);																	// calculate j1 for the given k
	i = (k-j_mod)/size_v;																			// calculate i for the given k

	//intE = Int_E(U,i); intE1 = Int_E1st(U,i); intE2 = Int_E2nd(U,i);
  
	if(intE[i]>0)																					// do this if the average direction of the field E over the space cell i is positive
	{
		j1r=j1+1;  j1l=j1;																			// set j1r to the value of j1+1 and j1l to the value of j1 (as here the the average flow of the field is from left to right so that gh^- must be used at the cell edges, as information flows against the field)
		kkr=i*size_v + (j1r*Nv*Nv + j2*Nv + j3);													// calculate the value of kkr for this value of j1r
		kkl=k; 																						// set kkl to k (since j1l = j1)
		if(j1r<Nv)ur = -U[kkr*6+2];																	// if j1r is not Nv (so that this cell is not receiving information from the right boundary), set ur to the negative of the coefficient of the basis function with shape 2 which is non-zero in the cell with global index kkr (which corresponds to the evaluation of gh^- at the right boundary and -ve since the field is negative here?) - note that if the cell was receiving information from the right boundary then gh^- = 0 here so ur is not needed
		ul = -U[kkl*6+2];																			// set ul to the negative of the coefficient of the basis function with shape 2 which is non-zero in the cell with global index kkl (which corresponds to the evaluation of gh^- at the right boundary and -ve since phi < 0 here)
	}
	else																							// do this if the average direction of the field E over the space cell i is non-positive
	{
		j1r=j1; j1l=j1-1;																			// set j1r to the value of j1 and j1l to the value of j1-1 (as here the the average flow of the field is from right to left so that gh^+ must be used at the cell edges, as information flows against the field)
		kkr=k;																						// set kkr to k (since j1r = j1)
		kkl=i*size_v + (j1l*Nv*Nv + j2*Nv + j3);													// calculate the value of kkl for this value of j1l
		ur = U[kkr*6+2];																			// set ur to the the coefficient of the basis function with shape 2 which is non-zero in the cell with global index kkr (which corresponds to the evaluation of gh^+ at the left boundary and +ve since phi > 0 here)
		if(j1l>-1)ul = U[kkl*6+2];																	// if j1l is not -1 (so that this cell is not receiving information from the left boundary), set ul to the coefficient of the basis function with shape 2 which is non-zero in the cell with global index kkl (which corresponds to the evaluation of gh^+ at the left boundary and +ve since phi < 0 here and being subtracted?) - note that if the cell was receiving information from the left boundary then gh^+ = 0 here so ul is not needed
	}

	if(l==0)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 0 (i.e. constant) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1) result = dv*dv*(U[kkr*6+0] + 0.5*ur + U[kkr*6+5]*5./12.- U[kkl*6+0] - 0.5*ul - U[kkl*6+5]*5./12.)*intE[i] + dv*dv*(U[kkr*6+1]-U[kkl*6+1])*intE1[i];		// this is the value at an interior cell
		else if(j1r<Nv)result =   dv*dv*(U[kkr*6+0] + 0.5*ur + U[kkr*6+5]*5./12.)*intE[i] + dv*dv*U[kkr*6+1]*intE1[i];																	// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		else if(j1l>-1)result = dv*dv*(- U[kkl*6+0] - 0.5*ul - U[kkl*6+5]*5./12.)*intE[i] - dv*dv*U[kkl*6+1]*intE1[i];																	// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
	}
	if(l==1)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 1 (i.e. linear in x) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1) result = dv*dv*( (U[kkr*6+0] + 0.5*ur + U[kkr*6+5]*5./12. - U[kkl*6+0] - 0.5*ul - U[kkl*6+5]*5./12.)*intE1[i] + (U[kkr*6+1] - U[kkl*6+1])*intE2[i] );		// this is the value at an interior cell
		else if(j1r<Nv)result=dv*dv*( (U[kkr*6+0] + 0.5*ur + U[kkr*6+5]*5./12.)*intE1[i] + U[kkr*6+1]*intE2[i] );																		// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		else if(j1l>-1)result=dv*dv*( (- U[kkl*6+0] - 0.5*ul - U[kkl*6+5]*5./12.)*intE1[i] - U[kkl*6+1]*intE2[i] );																		// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
	}
	if(l==2)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 2 (i.e. linear in v1) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1)result = 0.5*(dv*dv*(U[kkr*6+0] + 0.5*ur + U[kkr*6+5]*5./12.+ U[kkl*6+0] + 0.5*ul + U[kkl*6+5]*5./12.)*intE[i] + dv*dv*(U[kkr*6+1]+U[kkl*6+1])*intE1[i]);	// this is the value at an interior cell
		else if(j1r<Nv)result = 0.5*(dv*dv*(U[kkr*6+0] + 0.5*ur + U[kkr*6+5]*5./12.)*intE[i] + dv*dv*U[kkr*6+1]*intE1[i]);																// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		else if(j1l>-1)result = 0.5*(dv*dv*(U[kkl*6+0] + 0.5*ul + U[kkl*6+5]*5./12.)*intE[i] + dv*dv*U[kkl*6+1]*intE1[i]);    															// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
	}
	if(l==3)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 3 (i.e. linear in v2) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1)result = (U[kkr*6+3]-U[kkl*6+3])*intE[i]*dv*dv/12.;						// this is the value at an interior cell
		else if(j1r<Nv)result = U[kkr*6+3]*intE[i]*dv*dv/12.;										// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		else if(j1l>-1)result = -U[kkl*6+3]*intE[i]*dv*dv/12.;										// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
	}
	if(l==4)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 4 (i.e. linear in v3) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1)result = (U[kkr*6+4]-U[kkl*6+4])*intE[i]*dv*dv/12.;						// this is the value at an interior cell
		else if(j1r<Nv)result=U[kkr*6+4]*intE[i]*dv*dv/12.;											// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		else if(j1l>-1)result=-U[kkl*6+4]*intE[i]*dv*dv/12.;										// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
	}
	if(l==5)																						// calculate \int_i E*f*phi dx at interface v1==v_j+1/2 - \int_i E*f*phi dx at interface v1==v_j-1/2 for the basis function with shape 5 (i.e. the modulus of v) which is non-zero in the cell with global index k
	{
		if(j1r<Nv && j1l>-1)result = dv*dv*( ((U[kkr*6+0] + 0.5*ur - U[kkl*6+0] - 0.5*ul)*5./12. + (U[kkr*6+5]- U[kkl*6+5])*133./720.)*intE[i] + (U[kkr*6+1] - U[kkl*6+1])*intE1[i]*5./12. ); //BUG: coefficient of U[k][5] was 11/48 insteadof 133/720		// this is the value at an interior cell
		else if(j1r<Nv)result= dv*dv*( ((U[kkr*6+0] + 0.5*ur)*5./12. + U[kkr*6+5]*133./720.)*intE[i] + U[kkr*6+1]*intE1[i]*5./12. );													// this is the value at the cell at the left boundary in v1 (so that the integral over the left edge is zero)
		else if(j1l>-1)result=-dv*dv*( ((U[kkl*6+0] + 0.5*ul)*5./12. + U[kkl*6+5]*133./720.)*intE[i] + U[kkl*6+1]*intE1[i]*5./12. );													// this is the value at the cell at the right boundary in v1 (so that the integral over the right edge is zero)
	}

  	return result;
}

double f_marg(double *U, int i, int j1, double x, double v1)
{
	int j2, j3, k0, j2N, k;																		// declare j2, j3 (the indices of the velocity cell in the v2 & v3 directions), k0 (to store the value of i*Nv^3 + j1_Nv^2), j2N (to store the value of j2*N) & k (the index of the cell in U)
	double x_dif, v1_dif, retn;																	// declare x_dif (to store x - x_i), v1_dif (to store v1 - v_j1) & retn (the value of the marginal evaluated at the given x & v1 to be returned at the end

	x_dif = x -  Gridx((double)i);																		// set x_dif to x - x_i
	v1_dif = v1 - Gridv((double)j1);																	// set v1_dif to v1 - v_j1
	k0 = i*Nv*Nv*Nv + j1*Nv*Nv;																	// set k0 to i*Nv^3 + j1*Nv^2
	retn = 0;																					// initialise fM_val at 0, since the value is calculated through a sum

	//#pragma omp parallel for schedule(dynamic) private(j2,j2N,j3,k)  shared(U,Nv,k0,dx,dv,x_dif,v1_dif) reduction(+:retn)
	for(j2=0; j2<Nv; j2++)
	{
		j2N = j2*Nv;																			// set j2N to j2*Nv
		for(j3=0; j3<Nv; j3++)
		{
			k = k0 + j2N + j3;																	// set k to i*Nv^3 + j1*Nv^2 + j2*Nv + j3
			retn += dv*dv*U[6*k+0] + dv*dv*U[6*k+1]*x_dif/dx + dv*U[6*k+2]*v1_dif
							+ U[6*k+5]*(v1_dif*v1_dif +dv*dv/6);								// add dv*dv*U[6*k+0] + dv*dv*U[6*k+1]*x_dif/dx + dv*U[6*k+2]*v1_dif + U[6*k+5]*(v1_dif*v1_dif +dv*dv/6) for the given j2 & j3 in the sum for retn
		}
	}
	return retn;																				// return the value of the marginal evaluated at x & v1
}

void PrintMarginalLoc(FILE *margfile)															// function to print the values of x & v1 which the marginal will be evaluated at in the first two rows of the file with tag margfile (subsequent rows will be the values of the marginal at given timesteps)
{
	int i, j1, np, nx, nv;																		// declare i (the index of the space cell),  j1 (the index of the velocity cell in the v1 direction), np (the number of points to evaluate in a given space/velocity cell), nx (a counter for the points in the space cell) & nv (a counter for the points in the velocity cell)
	double x_0, x_val, v1_0, v1_val, ddx, ddv;													// declare x_0 (the x value at the left edge of a given cell), x_val (the x value to be evaluated at), declare v1_0 (the v1 value at the left edge of a given cell), v1_val (the v1 value to be evaluated at), ddx (the space between x values) & ddv (the space between v1 values)

	np = 4;																						// set np to 4
	ddx = dx/np;																				// set ddx to the space cell width divided by np
	ddv = dv/np;																				// set ddv to the velocity cell width divided by np
	for(i=0; i<Nx; i++)
	{
		x_0 = Gridx((double)i - 0.5);															// set x_0 to the value of x at the left edge of the i-th space cell
		for (nx=0; nx<np; nx++)
		{
			x_val = x_0 + nx*ddx;																// set x_val to x_0 plus nx increments of width ddx
			for(j1=0; j1<Nv; j1++)
			{
				for (nv=0; nv<np; nv++)
				{
					fprintf(margfile, "%11.8g  ", x_val);												// in the file tagged as fmarg, print the x coordinate
				}
			}
		}
	}
	fprintf(margfile, "\n");																	// print a new line in the file tagged as fmarg
	for(i=0; i<Nx; i++)
	{
		for (nx=0; nx<np; nx++)
		{
			for(j1=0; j1<Nv; j1++)
			{
				v1_0 = Gridv((double)j1 - 0.5);															// set v1_0 to the value of v1 at the left edge of the j1-th velocity cell in the v1 direction
				for (nv=0; nv<np; nv++)
				{
					v1_val = v1_0 + nv*ddv;																// set v1_val to v1_0 plus nv increments of width ddv
					fprintf(margfile, "%11.8g  ", v1_val);												// in the file tagged as fmarg, print the v1 coordinate
				}
			}
		}
	}
	fprintf(margfile, "\n");																	// print a new line in the file tagged as fmarg
}

void PrintMarginal(double *U, FILE *margfile)													// function to print the values of the marginal in the file tagged as margfile at the given timestep
{
	int i, j1, np, nx, nv;																		// declare i (the index of the space cell),  j1 (the index of the velocity cell in the v1 direction), np (the number of points to evaluate in a given space/velocity cell), nx (a counter for the points in the space cell) & nv (a counter for the points in the velocity cell)
	double x_0, x_val, v1_0, v1_val, fM_val, ddx, ddv;											// declare x_0 (the x value at the left edge of a given cell), x_val (the x value to be evaluated at), declare v1_0 (the v1 value at the left edge of a given cell), v1_val (the v1 value to be evaluated at), fM_val (the value of the marginal evaluated at (x_val, v1_val), ddx (the space between x values) & ddv (the space between v1 values)

	np = 4;																						// set np to 4
	ddx = dx/np;																				// set ddx to the space cell width divided by np
	ddv = dv/np;																				// set ddv to the velocity cell width divided by np
	for(i=0; i<Nx; i++)
	{
		x_0 = Gridx((double)i - 0.5);															// set x_0 to the value of x at the left edge of the i-th space cell
		for (nx=0; nx<np; nx++)
		{
			x_val = x_0 + nx*ddx;																// set x_val to x_0 plus nx increments of width ddx
			for(j1=0; j1<Nv; j1++)
			{
				v1_0 = Gridv((double)j1 - 0.5);													// set v1_0 to the value of v1 at the left edge of the j1-th velocity cell in the v1 direction
				for (nv=0; nv<np; nv++)
				{
					v1_val = v1_0 + nv*ddv;														// set v1_val to v1_0 plus nv increments of width ddv

					fM_val = f_marg(U, i, j1, x_val, v1_val);									// calculate the value of the marginal, evaluated at x_val & v1_val by using the function in the space cell i and velocity cell j1 in the v1 direction
					fprintf(margfile, "%11.8g  ", fM_val);										// in the file tagged as fmarg, print the value of the marginal f_M(t, x, v1)
				}
			}
		}
	}
	fprintf(margfile, "\n");																	// print a new line in the file tagged as fmarg
}

void ExportRhoQuadVals(double *U)																// function to export the values of the density rho suitable for Gaussian quadrature (should be done at equilibrium)
{
	int i, nx;																					// declare i (the index of the space cell) & nx (a counter for the space cell)
	double x_0, x_val;																	// declare x_0 (to store the x coordinate in the middle of the current cell), x_val (to store the x value to be evaluated at) & rho_val (to store the value of the density rho evaluated at the current x value)
	double *rho_vals;																			// declare a pointer to rho_vals (where the values of the density rho at equilibrium will be stored)
	rho_vals = (double*)malloc(5*Nx*sizeof(double));											// allocate enough space at the pointer rho_vals for 5*Nx many double numbers

	char flag[100], buffer_rhoeq[100];															// declare flag (a flag at the end of the file name, of 100 characters) & buffer_rhoeq (to store the name of the file where rho is store)
	sprintf(flag,"4Hump");																		// store the string "4Hump" in flag
	sprintf(buffer_rhoeq,"Data/RhoEquiVals_nu%gA%gk%gNx%dLx%gNv%dLv%gSpectralN%ddt%gnT%d_%s.dc",
					nu, A_amp, k_wave, Nx, Lx, Nv, Lv, N, dt, nT, flag);						// create a .dc file name, located in the directory Data, whose name is RhoEquiVals_ followed by the values of nu, A_amp, k_wave, Nx, Lx, Nv, Lv, N, dt, nT and the contents of flag and store it in buffer_rhoeq
	FILE *rhoeqvals;																			// declare a pointer to a file called rhoeqvals
	rhoeqvals = fopen(buffer_rhoeq,"w");														// set rhoeqvals to be a file with the name stored in buffer_rhoeq and set the file access mode of rhoeqvals to w (which creates an empty file and allows it to be written to)

	for(i=0;i<Nx;i++)																			// loop through the space cells
	{
		x_0 = Gridx((double)i);																	// set x_0 to the value of x at the center of the ith space cell

		for(nx=0;nx<5;nx++)																		// loop through the five quadrature points in the x direction of the cell
		{
			x_val = x_0 + 0.5*vt[nx]*dx;														// set x_val to the nx-th quadrature point in the cell
			rho_vals[i*5+nx] = rho_x(x_val, U, i);												// calculate the value of rho, evaluated at x_val by using the function in the space cell
		}
	}
	fwrite(rho_vals,sizeof(double),5*Nx,rhoeqvals);												// write the values of the density, stored in rhovals, which is 5*Nx entries, each of the size of a double datatype, in the file tagged as rhoeqvals
	fclose(rhoeqvals);																			// close the file rhoeqvals
}

void ComputeEquiVals(double *f)																	// function to compute the equilibrium solution, suitable for Gaussian quadrature and save it in the function f
{

	int i, j1, j2, j3, iNNN, j1NN, j2N, k, k_equi, nx, nv1, nv2, nv3;							// declare i (the index of the space cell), j1, j2, j3 (the indices of the velocity cell), iNNN (to store i*Nv^3), j1NN (to store j1*Nv^2), j2N (to store j2*Nv), k (the location of the given cell in U), k_equi (the location of the current quadrature point in f_equi), nx (a counter for the space cell) & nv1, nv2, nv3 (counters for the velocity cell)
	double x_0, v1_0, v2_0, v3_0, x_val, v1_val, v2_val, v3_val;								// declare x_0, v1_0, v2_0, v3_0 (to store the coordinates in the middle of the current cell), x_val, v1_val, v2_val & v3_val (to store the x & v values to be evaluated at)
	double *rho_vals;																			// declare a pointer to rho_vals (where the values of the density rho at equilibrium will be stored)
	rho_vals = (double*)malloc(5*Nx*sizeof(double));											// allocate enough space at the pointer rho_vals for 5*Nx many double numbers
	FILE *rhoeqvals;																			// declare a pointer to a file called rhoeqvals
	rhoeqvals = fopen("Data/RhoEquiVals_nu0.05A0.2k0.5Nx24Lx12.5664Nv24Lv5.25SpectralN16dt0.01nT0_4Hump.dc","r");														// set rhoeqvals to be a file with the name stored in buffer_rhoeq and setthe file access more of rhoeqvals to r (which allows the file to be read from)
	fread(rho_vals, sizeof(double), 5*Nx, rhoeqvals);											// read from the file rhoeqvals, which contains 5*Nx many entries of the size of a double number and store it rho_vals
	fclose(rhoeqvals);																			// close the file rhoeqvals
	for(i=0;i<5*Nx;i++)
	{
		printf("%g ", rho_vals[i]);
	}
	#pragma omp parallel for private(k,k_equi,i,j1,j2,j3,iNNN,j1NN,j2N,nx,nv1,nv2,nv3,x_0,v1_0,v2_0,v3_0,x_val,v1_val,v2_val,v3_val) shared(size_v,Nv,dv,dx,vt,wt,rho_vals,f)
	for(i=0;i<Nx;i++)																												// loop through the space cells
	{
		iNNN = i*size_v;																		// set iNNN to i*Nv^3
		x_0 = Gridx((double)i);																	// set x_0 to the value of x at the center of the ith space cell

		for(j1=0;j1<Nv;j1++)																	// loop through the velocity cells in the v1 direction
		{
			j1NN = j1*Nv*Nv;																	// set j1NN to j1*Nv^2
			v1_0 = Gridv((double)j1);															// set v1_0 to the value of v1 at the center of the j1th velocity cell in that direction
			for(j2=0;j2<Nv;j2++)																// loop through the velocity cells in the v2 direction
			{
				j2N = j2*Nv;																	// set j2N to j2*Nv
				v2_0 = Gridv((double)j2);														// set v2_0 to the value of v2 at the center of the j2th velocity cell in that direction
				for(j3=0;j3<Nv;j3++)															// loop through the velocity cells in the v3 direction
				{
					k = iNNN + j1NN + j2N + j3;													// set k to i*Nv^3 + j1*Nv^2 + j2*Nv + j3
					v3_0 = Gridv((double)j3);													// set v3_0 to the value of v3 at the center of the j3th velocity cell in that direction
					for(nx=0;nx<5;nx++)															// loop through the five quadrature points in the x direction of the cell
					{
						x_val = x_0 + 0.5*vt[nx]*dx;											// set x_val to the nx-th quadrature point in the cell
						for(nv1=0;nv1<5;nv1++)													// loop through the five quadrature points in the v1 direction of the cell
						{
							v1_val = v1_0 + 0.5*vt[nv1]*dv;										// set v1_val to the nv1-th quadrature point in the cell
							for(nv2=0;nv2<5;nv2++)												// loop through the five quadrature points in the v2 direction of the cell
							{
								v2_val = v2_0 + 0.5*vt[nv2]*dv;									// set v2_val to the nv2-th quadrature point in the cell
								for(nv3=0;nv3<5;nv3++)											// loop through the five quadrature points in the v3 direction of the cell
								{
									v3_val = v3_0 + 0.5*vt[nv3]*dv;								// set v3_val to the nv3-th quadrature point in the cell

									k_equi = 5*(5*(5*(5*k + nx) + nv1) + nv2) + nv3;			// set k_equi to the location that this value of f should be stored in f
									f[k_equi] = rho_vals[5*i+nx]*exp(-(v1_val*v1_val+v2_val*v2_val+v3_val*v3_val)/(2.1*2))/(sqrt(2.1*2*PI)*(2.1*2*PI));					// calculate the value of the equilbrium marginal, evaluated at x_val & v1_val by using the function in the space cell i and velocity cell (j1,j2,j3), namely rho_val(x)*exp(-(v1^2+v2^2+v3^2)/(2*T))/(2*T*PI)^(3/2)
								}
							}
						}
					}
				}
			}
		}
	}
}

void PrintEquiVals(double *U, FILE *margfile)
{
	int i, j1, np, nx, nv;																		// declare i (the index of the space cell),  j1 (the index of the velocity cell in the v1 direction), np (the number of points to evaluate in a given space/velocity cell), nx (a counter for the points in the space cell) & nv (a counter for the points in the velocity cell)
	double x_0, x_val, v1_0, v1_val, fM_val, rho_val, ddx, ddv;									// declare x_0 (the x value at the left edge of a given cell), x_val (the x value to be evaluated at), declare v1_0 (the v1 value at the left edge of a given cell), v1_val (the v1 value to be evaluated at), fM_val (the value of the marginal evaluated at (x_val, v1_val), rho_val (the value of the density evaluated at x_val), ddx (the space between x values) & ddv (the space between v1 values)

	np = 4;																						// set np to 4
	ddx = dx/np;																				// set ddx to the space cell width divided by np
	ddv = dv/np;																				// set ddv to the velocity cell width divided by np
	for(i=0; i<Nx; i++)
	{
		x_0 = Gridx((double)i - 0.5);															// set x_0 to the value of x at the left edge of the i-th space cell
		for (nx=0; nx<np; nx++)
		{
			x_val = x_0 + nx*ddx;																// set x_val to x_0 plus nx increments of width ddx
			for(j1=0; j1<Nv; j1++)
			{
				for (nv=0; nv<np; nv++)
				{
					fprintf(margfile, "%11.8g  ", x_val);												// in the file tagged as fmarg, print the x coordinate
				}
			}
		}
	}
	fprintf(margfile, "\n");																	// print a new line in the file tagged as fmarg
	for(i=0; i<Nx; i++)
	{
		for (nx=0; nx<np; nx++)
		{
			for(j1=0; j1<Nv; j1++)
			{
				v1_0 = Gridv((double)j1 - 0.5);															// set v1_0 to the value of v1 at the left edge of the j1-th velocity cell in the v1 direction
				for (nv=0; nv<np; nv++)
				{
					v1_val = v1_0 + nv*ddv;																// set v1_val to v1_0 plus nv increments of width ddv
					fprintf(margfile, "%11.8g  ", v1_val);												// in the file tagged as fmarg, print the v1 coordinate
				}
			}
		}
	}
	fprintf(margfile, "\n");																	// print a new line in the file tagged as fmarg
	for(i=0; i<Nx; i++)
	{
		x_0 = Gridx((double)i - 0.5);															// set x_0 to the value of x at the left edge of the i-th space cell
		for (nx=0; nx<np; nx++)
		{
			x_val = x_0 + nx*ddx;																// set x_val to x_0 plus nx increments of width ddx
			for(j1=0; j1<Nv; j1++)
			{
				v1_0 = Gridv((double)j1 - 0.5);													// set v1_0 to the value of v1 at the left edge of the j1-th velocity cell in the v1 direction
				for (nv=0; nv<np; nv++)
				{
					v1_val = v1_0 + nv*ddv;														// set v1_val to v1_0 plus nv increments of width ddv

					rho_val = rho_x(x_val, U, i);												// calculate the value of rho, evaluated at x_val by using the function in the space cell
					fM_val = rho_val*exp(-v1_val*v1_val/(2.1*2))/sqrt(2.1*2*PI);					// calculate the value of the marginal, evaluated at x_val & v1_val by using the function in the space cell i and velocity cell j1 in the v1 direction, namely rho_val*exp(-v1^2/(2*T))/sqrt(2*T*PI)
					fprintf(margfile, "%11.8g  ", fM_val);										// in the file tagged as fmarg, print the value of the marginal f_M(t, x, v1)
				}
			}
		}
	}
	fprintf(margfile, "\n");																	// print a new line in the file tagged as fmarg
}

double computeCellAvg(double *U, int i, int j1, int j2, int j3)																		// function to calculate the average value of the approximate function f (with DG coefficients in U) on the cell I_i x K_(j1,j2,j3), namely (1/cell_volume)*int_(I_i x K_(j1,j2,j3)) f dxdv = (1/(dx*dv^3))*int_(I_i x K_(j1,j2,j3)) f dxdv
{
	int k, nx, nv1, nv2, nv3;																										// declare k (the location of the given cell in U), nx (a counter for the space cell) & nv1, nv2, nv3 (counters for the velocity cell)
	double x_0, v1_0, v2_0, v3_0, x_val, v1_val, v2_val, v3_val, f_val, avg;														// declare x_0, v1_0, v2_0, v3_0 (to store the values in the middle of the current cell), x_val, v1_val, v2_val, v3_val (to store the x & v values to be evaluated at), f_val (to store the value of the function evaluated at the current x & v values) & avg (to store the average to be returned)
	k = Nv*(Nv*(i*Nv + j1) + j2) + j3;																								// set k to i*Nv^3 + j1*Nv^2 + j2*Nv + j3
	x_0 = Gridx((double)i);																											// set x_0 to the value of x at the center of the ith space cell
	v1_0 = Gridv((double)j1);																										// set v1_0 to the value of v1 at the center of the j1th velocity cell in that direction
	v2_0 = Gridv((double)j2);																										// set v2_0 to the value of v2 at the center of the j2th velocity cell in that direction
	v3_0 = Gridv((double)j3);																										// set v3_0 to the value of v3 at the center of the j3th velocity cell in that direction
	avg = 0;																														// set avg to 0 originally

	//#pragma omp parallel for private(nx,nv1,nv2,nv3,x_val,v1_val,v2_val,v3_val,f_val) shared(U,dv,dx,vt,wt,x_0,v1_0,v2_0,v3_0,k) reduction(+:avg)
	for(nx=0;nx<5;nx++)																												// loop through the five quadrature points in the x direction of the cell
	{
		x_val = x_0 + 0.5*vt[nx]*dx;																								// set x_val to the nx-th quadrature point in the cell
		for(nv1=0;nv1<5;nv1++)																										// loop through the five quadrature points in the v1 direction of the cell
		{
			v1_val = v1_0 + 0.5*vt[nv1]*dv;																							// set x_val to the nv1-th quadrature point in the cell
			for(nv2=0;nv2<5;nv2++)																									// loop through the five quadrature points in the v2 direction of the cell
			{
				v2_val = v2_0 + 0.5*vt[nv2]*dv;																						// set x_val to the nv2-th quadrature point in the cell
				for(nv3=0;nv3<5;nv3++)																								// loop through the five quadrature points in the v3 direction of the cell
				{
					v3_val = v3_0 + 0.5*vt[nv3]*dv;																					// set x_val to the nv3-th quadrature point in the cell
					f_val = U[k*6+0] + U[k*6+1]*(x_val-x_0)/dx + U[k*6+2]*(v1_val-v1_0)/dv +
							U[k*6+3]*(v2_val-v2_0)/dv + U[k*6+4]*(v3_val-v3_0)/dv +
							U[k*6+5]*(((v1_val-v1_0)/dv)*((v1_val-v1_0)/dv)+((v2_val-v2_0)/dv)*((v2_val-v2_0)/dv)
									+((v3_val-v3_0)/dv)*((v3_val-v3_0)/dv));														// set f_val to the evaluation of the approximation at (x_val,v1_val,v2_val,v3_val)
					avg += wt[nx]*wt[nv1]*wt[nv2]*wt[nv3]*f_val;																	// add f(x_val,v1_val,v2_val,v3_val) times the quadrature weights corresponding to the current x & v values to the current quadrature value
				}
			}
		}
	}
	avg = avg*0.5*0.5*0.5*0.5;																										// multiply the quadrature result by 0.5*0.5*0.5*0.5, since each quadrature is done over intervals of width dv (three times) or dx (once) instead of the standard interval [-1,1] (should be multiplied by dx*dv^3 also but then, to find the average, need to divide by dx*dv^3)
	return avg;																														// return the average value of f on the cell
}

void FindNegVals(double *U, int *NegVals, double *AvgVals)																							// function to find out the cells in which the approximation from U is negative on average and stores the cell locations in NegVals
{
	int i, j1, j2, j3, iNNN, j1NN, j2N, k, nx, nv1, nv2, nv3;																		// declare i (the index of the space cell), j1, j2, j3 (the indices of the velocity cell), iNNN (to store i*Nv^3), j1NN (to store j1*Nv^2), j2N (to store j2*Nv), k (the location of the given cell in U), nx (a counter for the space cell) & nv1, nv2, nv3 (counters for the velocity cell)
	double x_val, v1_val, v2_val, v3_val, f_avg;																					// declare x_val, v1_val, v2_val, v3_val (to store the x & v values to be evaluated at) & f_val (to store the value of the function evaluated at the current x & v values)
	#pragma omp parallel for private(i,j1,j2,j3,k,iNNN,j1NN,j2N,f_avg) shared(U,NegVals,AvgVals,size_v,Nv)
	for(i=0;i<Nx;i++)																												// loop through the space cells
	{
		iNNN = i*size_v;																											// set iNNN to i*Nv^3
		for(j1=0;j1<Nv;j1++)																										// loop through the velocity cells in the v1 direction
		{
			j1NN = j1*Nv*Nv;																										// set j1NN to j1*Nv^2
			for(j2=0;j2<Nv;j2++)																									// loop through the velocity cells in the v2 direction
			{
				j2N = j2*Nv;																										// set j2N to j2*Nv
				for(j3=0;j3<Nv;j3++)																								// loop through the velocity cells in the v3 direction
				{
					k = iNNN + j1NN + j2N + j3;																						// set k to i*Nv^3 + j1*Nv^2 + j2*Nv + j3
					NegVals[k] = 0;																									// set NegVals[k] to 0, which assumes at first that there will be no negative values in the kth cell
					f_avg = computeCellAvg(U,i,j1,j2,j3);																			// calculate the average value of the approximate solution in the current cell and set it to f_avg
					AvgVals[k] = f_avg;																								// store f_avg in AvgVals[k]
					if(f_avg < 0)																									// check if this value was negative
					{
						NegVals[k] = 1;																								// if so, set NegVals[k] to 1 to indicate that there was a negative value in this cell
					}
				}
			}
		}
	}
}

void FindNegVals_old(double *U, int *NegVals)																							// function to find out the cells in which the approximation from U turns negative and stores the cell locations in NegVals
{
	int i, j1, j2, j3, iNNN, j1NN, j2N, k, nx, nv1, nv2, nv3;																		// declare i (the index of the space cell), j1, j2, j3 (the indices of the velocity cell), iNNN (to store i*Nv^3), j1NN (to store j1*Nv^2), j2N (to store j2*Nv), k (the location of the given cell in U), nx (a counter for the space cell) & nv1, nv2, nv3 (counters for the velocity cell)
	double x_0, v1_0, v2_0, v3_0, x_val, v1_val, v2_val, v3_val, f_val;																// declare x_0, v1_0, v2_0, v3_0 (to store the coordinates in the middle of the current cell), x_val, v1_val, v2_val, v3_val (to store the x & v values to be evaluated at) & f_val (to store the value of the function evaluated at the current x & v values)
	for(i=0;i<Nx;i++)																												// loop through the space cells
	{
		iNNN = i*size_v;																											// set iNNN to i*Nv^3
		x_0 = Gridx((double)i);																												// set x_0 to the value of x at the center of the ith space cell

		for(j1=0;j1<Nv;j1++)																										// loop through the velocity cells in the v1 direction
		{
			j1NN = j1*Nv*Nv;																										// set j1NN to j1*Nv^2
			v1_0 = Gridv((double)j1);																										// set v1_0 to the value of v1 at the center of the j1th velocity cell in that direction
			for(j2=0;j2<Nv;j2++)																									// loop through the velocity cells in the v2 direction
			{
				j2N = j2*Nv;																										// set j2N to j2*Nv
				v2_0 = Gridv((double)j2);																									// set v2_0 to the value of v2 at the center of the j2th velocity cell in that direction
				for(j3=0;j3<Nv;j3++)																								// loop through the velocity cells in the v3 direction
				{
					k = iNNN + j1NN + j2N + j3;																						// set k to i*Nv^3 + j1*Nv^2 + j2*Nv + j3
					v3_0 = Gridv((double)j3);																								// set v3_0 to the value of v3 at the center of the j3th velocity cell in that direction
					NegVals[k] = 0;																									// set NegVals[k] to 0, which assumes at first that there will be no negative values in the kth cell
					for(nx=0;nx<5;nx++)																								// loop through the five quadrature points in the x direction of the cell
					{
						x_val = x_0 + 0.5*vt[nx]*dx;																					// set x_val to the nx-th quadrature point in the cell
						for(nv1=0;nv1<5;nv1++)																						// loop through the five quadrature points in the v1 direction of the cell
						{
							v1_val = v1_0 + 0.5*vt[nv1]*dv;																			// set x_val to the nv1-th quadrature point in the cell
							for(nv2=0;nv2<5;nv2++)																					// loop through the five quadrature points in the v2 direction of the cell
							{
								v2_val = v2_0 + 0.5*vt[nv2]*dv;																		// set x_val to the nv2-th quadrature point in the cell
								for(nv3=0;nv3<5;nv3++)																				// loop through the five quadrature points in the v3 direction of the cell
								{
									v3_val = v3_0 + 0.5*vt[nv3]*dv;																	// set x_val to the nv3-th quadrature point in the cell
									f_val = U[k*6+0] + U[k*6+1]*(x_val-x_0)/dx + U[k*6+2]*(v1_val-v1_0)/dv +
											U[k*6+3]*(v2_val-v2_0)/dv + U[k*6+4]*(v3_val-v3_0)/dv +
											U[k*6+5]*(((v1_val-v1_0)/dv)*((v1_val-v1_0)/dv)+((v2_val-v2_0)/dv)*((v2_val-v2_0)/dv)
													+((v3_val-v3_0)/dv)*((v3_val-v3_0)/dv));										// set f_val to the evaluation of the approximation at (x_val,v1_val,v2_val,v3_val)
									if(f_val < 0)																					// check if this value was negative
									{
										NegVals[k] = 1;																				// if so, set NegVals[k] to 1 to indicate that there was a negative value in this cell
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void CheckNegVals(double *U, int *NegVals, double *AvgVals)																							// function to find out the cells in which the approximation from U turns negative and stores the cell locations in NegVals
{
	int i, j1, j2, j3, iNNN, j1NN, j2N, k;																							// declare i (the index of the space cell), j1, j2, j3 (the indices of the velocity cell), iNNN (to store i*Nv^3), j1NN (to store j1*Nv^2), j2N (to store j2*Nv), k (the location of the given cell in U), nx (a counter for the space cell) & nv1, nv2, nv3 (counters for the velocity cell)
	double f_avg;																													// declare x_0, v1_0, v2_0, v3_0 (to store the coordinates in the middle of the current cell), x_val, v1_val, v2_val, v3_val (to store the x & v values to be evaluated at) & f_val (to store the value of the function evaluated at the current x & v values)
	for(i=0;i<Nx;i++)
	{
		iNNN = i*size_v;																											// set iNNN to i*Nv^3
		for(j1=0;j1<Nv;j1++)																										// loop through the velocity cells in the v1 direction
		{
			j1NN = j1*Nv*Nv;																										// set j1NN to j1*Nv^2
			for(j2=0;j2<Nv;j2++)																									// loop through the velocity cells in the v2 direction
			{
				j2N = j2*Nv;																										// set j2N to j2*Nv
				for(j3=0;j3<Nv;j3++)																								// loop through the velocity cells in the v3 direction
				{
					k = iNNN + j1NN + j2N + j3;																						// set k to i*Nv^3 + j1*Nv^2 + j2*Nv + j3
					printf("%d ~ (%d, %d, %d, %d) f_avg = %g \n", NegVals[k], i, j1, j2, j3, AvgVals[k]);
					printf(" \n");
				}
			}
		}
	}
}

#ifdef MPI
void computeH(double *H, double *U)// H_k(i,j)(f, E, phi_l)  
{
  int k, l; // k=i*Nv^3 + (j1*Nv*Nv + j2*Nv + j3)
  double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tp0, tp5;
 
  #pragma omp parallel for schedule(dynamic) private(tp0, tp5, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6) shared(H, U)
  for(k=chunksize_dg*myrank_mpi;k<chunksize_dg*(myrank_mpi+1);k++){        
	tp0 = I1(U,k,0) - I2(U,k,0) - I3(U,k,0) + I5(U,k,0);
	tp5 = I1(U,k,5) - I2(U,k,5) - I3(U,k,5) + I5(U,k,5);
	H[(k%chunksize_dg)*6+0] = (19*tp0/4. - 15*tp5)/(dx*dv*dv*dv);
	H[(k%chunksize_dg)*6+5] = (60*tp5 - 15*tp0)/(dx*dv*dv*dv);	
    for(l=1;l<5;l++){
	  tmp1=I1(U,k,l); tmp2=I2(U,k,l); tmp3=I3(U,k,l);  tmp5=I5(U,k,l); 
	  H[(k%chunksize_dg)*6+l] = (tmp1-tmp2-tmp3+tmp5)*12./(dx*dv*dv*dv);				  
	}	
  }  
}

void RK3(double *U) // RK3 for f_t = H(f)
{
  int i, k, l, k_local;
  double tp0, tp1, tp2, tp3, tp4, tp5, H[6];//, tp0, tp5, tmp1, tmp2, tmp3, tmp5;
 
  MPI_Status status;
  
  ce = computePhi_x_0(U); 

  #pragma omp parallel for private(i) shared(U,cp, intE, intE1)
  for(i=0;i<Nx;i++){
    cp[i] = computeC_rho(U,i); intE[i] = Int_E(U,i); intE1[i] = Int_E1st(U,i); 
  }
  
  #pragma omp parallel for private(i) shared(intE2)
  for(i=0;i<Nx;i++){
    intE2[i] = Int_E2nd(U,i); // BUG: Int_E2nd() require knowldege of cp 
  }
  
  #pragma omp parallel for schedule(dynamic)  private(H,k, k_local, l, tp0, tp1, tp2, tp3, tp4, tp5) shared(U, Utmp)
  for(k=chunksize_dg*myrank_mpi;k<chunksize_dg*(myrank_mpi+1);k++){ 
    k_local = k%chunksize_dg;
    
    tp0=I1(U,k,0)-I2(U,k,0)-I3(U,k,0)+I5(U,k,0);
    tp1=I1(U,k,1)-I2(U,k,1)-I3(U,k,1)+I5(U,k,1);
    tp2=I1(U,k,2)-I2(U,k,2)-I3(U,k,2)+I5(U,k,2);
    tp3=I1(U,k,3)-I2(U,k,3)-I3(U,k,3)+I5(U,k,3);
    tp4=I1(U,k,4)-I2(U,k,4)-I3(U,k,4)+I5(U,k,4);
    tp5=I1(U,k,5)-I2(U,k,5)-I3(U,k,5)+I5(U,k,5);

    //H[k_local][0] = (19*tp[0]/4. - 15*tp[5])/dx/scalev;
    //H[k_local][5] = (60*tp[5] - 15*tp[0])/dx/scalev;	
    H[0] = (19*tp0/4. - 15*tp5)/dx/scalev;
    H[5] = (60*tp5 - 15*tp0)/dx/scalev;	
    //for(l=1;l<5;l++)H[l] = tp[l]*12./dx/scalev;;//H[k_local][l] = tp[l]*12./dx/scalev;
    H[1] = tp1*12./dx/scalev; H[2] = tp2*12./dx/scalev; H[3] = tp3*12./dx/scalev; H[4] = tp4*12./dx/scalev;
    
    for(l=0;l<6;l++) Utmp[k_local*6+l] = U[k*6+l] + dt*H[l];	
  }    
  if(myrank_mpi == 0) {
    //dump the weights we've computed into U1
    for(k=0;k<chunksize_dg;k++) {
	for(l=0;l<6;l++) U1[k*6+l] = Utmp[k*6+l];	
    } 
    //receive from all other processes
    for(i=1;i<nprocs_mpi;i++) {
	    MPI_Recv(output_buffer_vp, chunksize_dg*6, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status); // receive weight from other processes consecutively, rank i=1..numNodes-1, ensuring the weights are stored in the file consecutively !
	    for(k=0;k<chunksize_dg;k++) {
			for(l=0;l<6;l++)U1[(k + i*chunksize_dg)*6+l] = output_buffer_vp[k*6+l];
		}			
    }
  }
  else  MPI_Send(Utmp, chunksize_dg*6, MPI_DOUBLE, 0, myrank_mpi, MPI_COMM_WORLD);
  
  
  MPI_Bcast(U1, size*6, MPI_DOUBLE, 0, MPI_COMM_WORLD);    
  MPI_Barrier(MPI_COMM_WORLD);
  /////////////////// 1st step of RK3 done//////////////////////////////////////////////////////// 
    
  ce = computePhi_x_0(U1); 

  #pragma omp parallel for private(i) shared(U1,cp, intE, intE1)
  for(i=0;i<Nx;i++){
    cp[i] = computeC_rho(U1,i); intE[i] = Int_E(U1,i); intE1[i] = Int_E1st(U1,i); 
  }
  
  #pragma omp parallel for private(i) shared(U1,intE2)
  for(i=0;i<Nx;i++){
    intE2[i] = Int_E2nd(U1,i); // BUG: Int_E2nd() require knowldege of cp 
  }
  
  #pragma omp parallel for schedule(dynamic) private(H, k, k_local, l, tp0, tp1, tp2, tp3, tp4, tp5)  shared(U,Utmp)
  for(k=chunksize_dg*myrank_mpi;k<chunksize_dg*(myrank_mpi+1);k++){      
    k_local = k%chunksize_dg;
    
    tp0=I1(U1,k,0)-I2(U1,k,0)-I3(U1,k,0)+I5(U1,k,0);
    tp1=I1(U1,k,1)-I2(U1,k,1)-I3(U1,k,1)+I5(U1,k,1);
    tp2=I1(U1,k,2)-I2(U1,k,2)-I3(U1,k,2)+I5(U1,k,2);
    tp3=I1(U1,k,3)-I2(U1,k,3)-I3(U1,k,3)+I5(U1,k,3);
    tp4=I1(U1,k,4)-I2(U1,k,4)-I3(U1,k,4)+I5(U1,k,4);
    tp5=I1(U1,k,5)-I2(U1,k,5)-I3(U1,k,5)+I5(U1,k,5);

    //H[k_local][0] = (19*tp[0]/4. - 15*tp[5])/dx/scalev;
    //H[k_local][5] = (60*tp[5] - 15*tp[0])/dx/scalev;	
    H[0] = (19*tp0/4. - 15*tp5)/dx/scalev;
    H[5] = (60*tp5 - 15*tp0)/dx/scalev;	
    //for(l=1;l<5;l++)H[l] = tp[l]*12./dx/scalev;;//H[k_local][l] = tp[l]*12./dx/scalev;
    H[1] = tp1*12./dx/scalev; H[2] = tp2*12./dx/scalev; H[3] = tp3*12./dx/scalev; H[4] = tp4*12./dx/scalev;
    
    for(l=0;l<6;l++) Utmp[k_local*6+l] = 0.75*U[k*6+l] + 0.25*U1[k*6+l] + 0.25*dt*H[l];
  }    
  if(myrank_mpi == 0) {
    //dump the weights we've computed into U1
    for(k=0;k<chunksize_dg;k++) {
	      for(l=0;l<6;l++) U1[k*6+l] = Utmp[k*6+l];	
    } 
    //receive from all other processes
    for(i=1;i<nprocs_mpi;i++) {      
		MPI_Recv(output_buffer_vp, chunksize_dg*6, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status); // receive weight from other processes consecutively, rank i=1..numNodes-1, ensuring the weights are stored in the file consecutively !
		for(k=0;k<chunksize_dg;k++) {
			for(l=0;l<6;l++)U1[(k + i*chunksize_dg)*6+l] = output_buffer_vp[k*6+l];
        }
    }
  }
  else MPI_Send(Utmp, chunksize_dg*6, MPI_DOUBLE, 0, myrank_mpi, MPI_COMM_WORLD);
  
  
  MPI_Bcast(U1, size*+6, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  /////////////////// 2nd step of RK3 done//////////////////////////////////////////////////////// 
   
ce = computePhi_x_0(U1); 

  #pragma omp parallel for private(i) shared(U1,cp, intE, intE1)
  for(i=0;i<Nx;i++){
    cp[i] = computeC_rho(U1,i); intE[i] = Int_E(U1,i); intE1[i] = Int_E1st(U1,i); 
  }
  
  #pragma omp parallel for private(i) shared(U1,intE2)
  for(i=0;i<Nx;i++){
    intE2[i] = Int_E2nd(U1,i); // BUG: Int_E2nd() require knowldege of cp 
  }
  
  #pragma omp parallel for schedule(dynamic) private(H, k, k_local, l, tp0, tp1, tp2, tp3, tp4, tp5)  shared(U,Utmp)
  for(k=chunksize_dg*myrank_mpi;k<chunksize_dg*(myrank_mpi+1);k++){      
    k_local = k%chunksize_dg;
  
    tp0=I1(U1,k,0)-I2(U1,k,0)-I3(U1,k,0)+I5(U1,k,0);
    tp1=I1(U1,k,1)-I2(U1,k,1)-I3(U1,k,1)+I5(U1,k,1);
    tp2=I1(U1,k,2)-I2(U1,k,2)-I3(U1,k,2)+I5(U1,k,2);
    tp3=I1(U1,k,3)-I2(U1,k,3)-I3(U1,k,3)+I5(U1,k,3);
    tp4=I1(U1,k,4)-I2(U1,k,4)-I3(U1,k,4)+I5(U1,k,4);
    tp5=I1(U1,k,5)-I2(U1,k,5)-I3(U1,k,5)+I5(U1,k,5);

    H[0] = (19*tp0/4. - 15*tp5)/dx/scalev;
    H[5] = (60*tp5 - 15*tp0)/dx/scalev;	
    //for(l=1;l<5;l++)H[l] = tp[l]*12./dx/scalev;;//H[k_local][l] = tp[l]*12./dx/scalev;
    H[1] = tp1*12./dx/scalev; H[2] = tp2*12./dx/scalev; H[3] = tp3*12./dx/scalev; H[4] = tp4*12./dx/scalev;	

    for(l=0;l<6;l++) Utmp[k_local*6+l] = U[k*6+l]/3. + U1[k*6+l]*2./3. + dt*H[l]*2./3.;
  }    
  if(myrank_mpi == 0) {
    //dump the weights we've computed into U1
    for(k=0;k<chunksize_dg;k++) {
	      for(l=0;l<6;l++) U[k*6+l] = Utmp[k*6+l];	
    } 
    //receive from all other processes
    for(i=1;i<nprocs_mpi;i++) {      
	    MPI_Recv(output_buffer_vp, chunksize_dg*6, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status); // receive weight from other processes consecutively, rank i=1..numNodes-1, ensuring the weights are stored in the file consecutively !
	    for(k=0;k<chunksize_dg;k++) {
		  for(l=0;l<6;l++)U[(k + i*chunksize_dg)*6+l] = output_buffer_vp[k*6+l];
        } 
    }
  }
  else MPI_Send(Utmp, chunksize_dg*6, MPI_DOUBLE, 0, myrank_mpi, MPI_COMM_WORLD);
  
  MPI_Bcast(U, size*6, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); 
  /////////////////// 3rd step of RK3 done//////////////////////////////////////////////////////// 
}


#else
void computeH(double *U)// H_k(i,j)(f, E, phi_l)  
{
  int i, k, l; // k=i*Nv^3 + (j1*Nv*Nv + j2*Nv + j3)
  double tp0, tp1, tp2, tp3, tp4, tp5;
 

  ce = computePhi_x_0(U); 
  // #pragma omp barrier
  #pragma omp parallel for private(i) shared(U,cp, intE, intE1)
  for(i=0;i<Nx;i++){
    cp[i] = computeC_rho(U,i); intE[i] = Int_E(U,i); intE1[i] = Int_E1st(U,i); 
  }
  
  #pragma omp parallel for private(i) shared(U, intE2)
  for(i=0;i<Nx;i++){
    intE2[i] = Int_E2nd(U,i); // BUG: Int_E2nd() require knowldege of cp 
  }
  
  #pragma omp parallel for schedule(dynamic) private(tp0, tp1, tp2, tp3, tp4, tp5, k,l) shared(U, H)
  for(k=0;k<size;k++){
   // double tp[6];

  // #pragma omp critical
   // {
      /*for(l=0;l<6;l++){
	tp[l]=I1(U,k,l)-I2(U,k,l)-I3(U,k,l)+I5(U,k,l);		
      } */
   // }
    tp0=I1(U,k,0)-I2(U,k,0)-I3(U,k,0)+I5(U,k,0);
    tp1=I1(U,k,1)-I2(U,k,1)-I3(U,k,1)+I5(U,k,1);
    tp2=I1(U,k,2)-I2(U,k,2)-I3(U,k,2)+I5(U,k,2);
    tp3=I1(U,k,3)-I2(U,k,3)-I3(U,k,3)+I5(U,k,3);
    tp4=I1(U,k,4)-I2(U,k,4)-I3(U,k,4)+I5(U,k,4);
    tp5=I1(U,k,5)-I2(U,k,5)-I3(U,k,5)+I5(U,k,5);
     
    H[k*6+0] = (19*tp0/4. - 15*tp5)/dx/scalev;
    H[k*6+5] = (60*tp5 - 15*tp0)/dx/scalev;	
    //for(l=1;l<5;l++)H[k][l] = tp[l]*12./dx/scalev;

    H[k*6+1] = tp1*12./dx/scalev; H[k*6+2] = tp2*12./dx/scalev; H[k*6+3] = tp3*12./dx/scalev; H[k*6+4] = tp4*12./dx/scalev;
  }
  
}

void RK3(double *U) // RK3 for f_t = H(f)
{
  int k, l;
  //double **H = (double **)malloc(size*sizeof(double *));
  //for (l=0;l<size;l++) H[l] = (double*)malloc(6*sizeof(double));
  //double **U1 = (double **)malloc(size*sizeof(double *));
  //for (l=0;l<size;l++) U1[l] = (double*)malloc(6*sizeof(double));
  
  computeH(U);
   #pragma omp parallel for private(k,l) shared(U)
  for(k=0;k<size;k++){				  
	  for(l=0;l<6;l++) U1[k*6+l] = U[k*6+l] + dt*H[k*6+l];				  
  }
 
  computeH(U1);
   #pragma omp parallel for private(k,l) shared(U)
  for(k=0;k<size;k++){						  
	  for(l=0;l<6;l++) U1[k*6+l] = 0.75*U[k*6+l] + 0.25*U1[k*6+l] + 0.25*dt*H[k*6+l];		 
  }

  computeH(U1);
   #pragma omp parallel for private(k,l) shared(U)
  for(k=0;k<size;k++){						  
	for(l=0;l<6;l++) U[k*6+l] = U[k*6+l]/3. + U1[k*6+l]*2./3. + dt*H[k*6+l]*2./3.;			
  }
  //free(H); free(U1);
}
#endif
