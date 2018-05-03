/* This is the source file which contains the subroutines necessary for calculating the electric potential and the
 * corresponding field, as well as any integrals which involve the field for the sake of the collisionless advection
 * problem.
 *
 * Functions included: Gridv, Gridx, rho_x, rho, computePhi_x_0, computePhi, PrintPhiVals, computeC_rho, Int_Int_rho
 * Int_Int_rho1st, Int_E, Int_E1st, Int_E2nd, Int_fE, I1, I2, I3, I5, computeH, RK3
 *
 */

#include "FieldCalculations.h"																	// FieldCalculations.h is where the prototypes for the functions contained in this file are declared

#ifdef Doping																					// only compile this file if Doping is defined and so the more complicated field calculations must be used

double DopingProfile(int i)																		// function to return a step function doping profile, based on what cell a given x is in
{
	double ND;																					// declare DP (the value of the doping profile in the current space cell to be returned)
	if(i <= a_i || i > b_i)																		// if the space cell i is either in the lower or upper third of all cells then set the value of DP to NH
	{
		ND = NH;
	}
	else																						// if the space cell i is either in the middle third of all cells then set the value of DP to NL
	{
		ND = NL;
	}
	return ND;																					// return the value of DP
}

void PrintDoping()
{
	double ND;
	printf("Doping Profile: ");
	for(int i=0; i<Nx; i++)
	{
		ND = DopingProfile(i);
		printf("ND[%d] = %g, ", i, ND);
	}
	printf("\n");
}

double EpsilonValue(int i)																		// function to return a step function doping profile, based on what cell a given x is in
{
	double epsilon;																					// declare DP (the value of the doping profile in the current space cell to be returned)
	if(i <= b_i)																		// if the space cell i is either in the lower or upper third of all cells then set the value of DP to NH
	{
		epsilon = 0.1;
	}
	else																						// if the space cell i is either in the middle third of all cells then set the value of DP to NL
	{
		epsilon = 0.1;
	}
	return epsilon;																					// return the value of DP
}

void PrintEpsilon()
{
	double eps;
	printf("Permittivity Values: ");
	for(int i=0; i<Nx; i++)
	{
		eps = EpsilonValue(i);
		printf("eps[%d] = %g, ", i, eps);
	}
	printf("\n");
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

double computePhi_x_0(double *U) /* DIFFERENT FOR withND */																// compute the constant coefficient of x in phi, which is actually phi_x(0) (Calculate C_E in the paper -between eq. 52 & 53?), with doping profile given by ND above
{
	int i, j, k, m, q;
	double tmp=0.;
	double a_val = (a_i+1)*dx;
	double b_val = (b_i+1)*dx;
	double Phi_Lx = -1;																									// declare Phi_Lx (the Dirichlet BC, Phi(t, L_x) = Phi_Lx) and set its value
	double eps;
	eps = EpsilonValue(Nx-1);																// set eps to the value of epsilon at ix

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

	#ifdef Electrons
	return Phi_Lx/Lx + 0.5*NH*Lx/eps + (NL-NH)*(b_val-a_val)/eps - (0.5*(NL-NH)*(b_val*b_val - a_val*a_val) + tmp)/(Lx*eps);
	#endif	// Electrons
	#ifdef Ions
	return Phi_Lx/Lx - (0.5*NH*Lx/eps + (NL-NH)*(b_val-a_val)/eps - (0.5*(NL-NH)*(b_val*b_val - a_val*a_val) + tmp)/(Lx*eps));
	#endif	// Ions
}

double computePhi(double *U, double x, int ix)	/* DIFFERENT FOR withND */											// function to compute the potential Phi at a position x, contained in [x_(ix-1/2), x_(ix+1/2)]
{
	int i_out, i, j1, j2, j3, iNNN, j1NN, j2N, k;										// declare counters i_out (for the outer sum of i values), i, j1, j2, j3 for summing the contribution from cell I_i x K_(j1, j2, j3), iNNN (the value of i*Nv^3), j1NN (the value of j1*Nv^2), j2N (the value of j2*N) & k (the location in U of the cell I_i x K_(j1, j2, j3))
	double retn, sum1, sum3, sum4, x_diff, x_diff_mid, x_diff_sq, x_eval, C_E;			// declare retn (the value of Phi returned at the end), sum1 (the value of the first two sums), sum3 (the value of the third sum), sum4 (the value of the fourth sum), x_diff (the value of x - x_(ix-1/2)), x_diff_mid (the value of x - x_ix), x_diff_sq (the value of x_diff^2), x_eval (the value associated to the integral of (x - x_i)^2) & C_E (the value of the constant in the formula for phi)
	double ND, eps;																			// declare ND (the value of the doping profile at the given x)
	sum1 = 0;
	sum3 = 0;
	sum4 = 0;
	retn = 0;
	ND = DopingProfile(ix);																// set ND to the value of the doping profile at ix
	eps = EpsilonValue(ix);																// set eps to the value of epsilon at ix
	x_diff = x - Gridx(ix-0.5);
	x_diff_mid = x - Gridx(ix);
	x_diff_sq = x_diff*x_diff;
	x_eval = x_diff_mid*x_diff_mid*x_diff_mid/(6.*dx) - dx*x_diff_mid/8. - dx*dx/24.;

	for(i_out = 0; i_out < ix; i_out++)
	{
		sum1 += computeC_rho(U, i_out);

		iNNN = i_out*Nv*Nv*Nv;
		#pragma omp parallel for private(j1,j2,j3,j1NN,j2N,k) shared(Nv,iNNN,U,dx,scalev) reduction(+:sum1)
		for(j1 = 0; j1 < Nv; j1++)
		{
			j1NN = j1*Nv*Nv;
			for(j2 = 0; j2 < Nv; j2++)
			{
				j2N = j2*Nv;
				for(j3 = 0; j3 < Nv; j3++)
				{
					k = iNNN + j1NN + j2N + j3;
					sum1 += (U[6*k]/2 - U[6*k+1]/12. +  U[6*k+5]/8)*dx*scalev;
				}
			}
		}
	}
	sum1 = sum1*dx;//*dx;

	sum3 = computeC_rho(U, ix);

	sum3 = sum3*x_diff;
	iNNN = ix*Nv*Nv*Nv;
	#pragma omp parallel for private(j1,j2,j3,j1NN,j2N,k) shared(Nv,iNNN,U,x_diff_sq,x_eval) reduction(+:sum4)
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
	sum4 = sum4*scalev;

	C_E = computePhi_x_0(U);
	retn = sum1 + sum3 + sum4 - ND*x*x/2.;
	if(ix > a_i)																						// if x > a then there is an extra term to add
	{
		double a_val = (a_i+1)*dx;																		// declare a_val and set it to the value of x at the edge of the ai-th space cell
		retn -= (NH-NL)*a_val*(x - 0.5*a_val);															// add (NH-NL)a(x-a/2) to retn
	}
	if(ix > b_i)																						// if x > b then there is an extra term to add
	{
		double b_val = (b_i+1)*dx;																		// declare b_val and set it to the value of x at the edge of the bi-th space cell
		retn -= (NL-NH)*b_val*(x - 0.5*b_val);															// add (NL-NH)b(x-b/2) to retn
	}

	#ifdef Electrons
	retn = retn/eps + C_E*x;
	#endif	// Electrons
	#ifdef Ions
	retn = -retn/eps + C_E*x;
	#endif	// Ions
	return retn;																						// return the value of phi at x
}

double computeE(double *U, double x, int ix)	/* DIFFERENT FOR withND */								// function to compute the field E at a position x, contained in [x_(ix-1/2), x_(ix+1/2)]
{
	int iNNN, j, k;																						// declare j (a counter for summing the contribution from the velocity cells), iNNN (the value of i*Nv^3) & k (the location in U of the cell I_ix x K_(j1, j2, j3))
	double retn, sum, x_diff, x_diff_mid, x_eval, C_E;													// declare retn (the value of E returned at the end), sum (the value of the sum to calculate the integral of rho), x_diff (the value of x - x_(ix-1/2)), x_diff_mid (the value of x - x_ix), x_eval (the value associated to the integral of (x - x_i)^2) & C_E (the value of the constant in the formula for E)
	double ND, eps;																							// declare ND (the value of the doping profile at the given x)
	sum = 0;																							// initialise the sum at 0
	ND = DopingProfile(ix);																				// set ND to the value of the doping profile at ix
	eps = EpsilonValue(ix);																				// set eps to the value of epsilon at ix
	C_E = computePhi_x_0(U);
	x_diff = x - Gridx(ix-0.5);
	x_diff_mid = x - Gridx(ix);
	x_eval = x_diff_mid*x_diff_mid/(2.*dx) - dx/8.;

	iNNN = ix*size_v;
	#pragma omp parallel for private(j,k) shared(size_v,iNNN,U,x_diff,x_eval) reduction(+:sum)
	for(j=0;j<size_v;j++)
	{
		k = iNNN + j;
		sum += U[6*k+0]*x_diff + U[6*k+1]*x_eval + U[6*k+5]*x_diff/4.;
	}
	sum = sum*scalev;
	sum += computeC_rho(U, ix);

	retn = ND*x - sum;
	if(ix > a_i)																						// if x > a then there is an extra term to add
	{
		double a_val = (a_i+1)*dx;																		// declare a_val and set it to the value of x at the edge of the ai-th space cell
		retn += (NH-NL)*a_val;																			// add (NH-NL)a to retn
	}
	if(ix > b_i)																						// if x > b then there is an extra term to add
	{
		double b_val = (b_i+1)*dx;																		// declare b_val and set it to the value of x at the edge of the bi-th space cell
		retn += (NL-NH)*b_val;																			// add (NL-NH)b to retn
	}

	#ifdef Electrons
	retn = retn/eps - C_E;
	#endif	// Electrons
	#ifdef Ions
	retn = -retn/eps - C_E;
	#endif	// Ions
	return retn;																						// return the value of phi at x

}

void PrintFieldLoc(FILE *phifile, FILE *Efile)							// function to print the values of x & v1 which the marginal will be evaluated at in the first two rows of the file with tag margfile (subsequent rows will be the values of the marginal at given timesteps)
{
	double x_0, x_val, ddx;												// declare x_0 (the x value at the left edge of a given cell), x_val (the x value to be evaluated at) & ddx (the space between x values)
	int np = 4;															// declare np (the number of points to evaluate in a given space cell) and set its value

	ddx = dx/np;														// set ddx to the space cell width divided by np
	for(int i=0; i<Nx; i++)
	{
		x_0 = Gridx((double)i - 0.5);									// set x_0 to the value of x at the left edge of the i-th space cell
		for(int nx=0; nx<np; nx++)
		{
			x_val = x_0 + nx*ddx;										// set x_val to x_0 plus nx increments of width ddx
			fprintf(phifile, "%11.8g  ", x_val);						// in the file tagged as phifile, print the x coordinate
			fprintf(Efile, "%11.8g  ", x_val);							// in the file tagged as Efile, print the x coordinate
		}
	}
	fprintf(phifile, "\n");												// print a new line in the file tagged as phifile
	fprintf(Efile, "\n");												// print a new line in the file tagged as Efile
}

void PrintFieldData(double* U_vals, FILE *phifile, FILE *Efile)							// function to print the values of the potential and the field in the x1 & x2 directions in the file tagged as phifile, Ex1file & Ex2file, respectively, at the given timestep
{
	double x_0, x_val, phi_val, E_val, ddx;								// declare x_0 (the x value at the left edge of a given cell), x_val (the x value to be evaluated at), phi_val (the value of phi evaluated at x_val), E_val (the value of E evaluated at x_val) & ddx (the space between x values)
	int np = 4;															// declare np (the number of points to evaluate in a given space cell) and set its value

	ddx = dx/np;														// set ddx to the space cell width divided by np
	for(int i=0; i<Nx; i++)
	{
		x_0 = Gridx((double)i - 0.5);									// set x_0 to the value of x at the left edge of the i-th space cell
		for (int nx=0; nx<np; nx++)
		{
			x_val = x_0 + nx*ddx;										// set x_val to x_0 plus nx increments of width ddx

			phi_val = computePhi(U_vals, x_val, i);							// calculate the value of phi, evaluated at x_val by using the function in the space cell i
			E_val = computeE(U_vals, x_val, i);								// calculate the value of E, evaluated at x_val by using the function in the space cell i
			fprintf(phifile, "%11.8g ", phi_val);						// in the file tagged as phifile, print the value of the potential phi(t, x_val)
			fprintf(Efile, "%11.8g ", E_val);							// in the file tagged as Efile, print the value of the field E(t, x_val)
		}
	}
	fprintf(phifile, "\n");												// in the file tagged as phifile, print a new line
	fprintf(Efile, "\n");												// in the file tagged as Efile, print a new line
}

void PrintPhiVals(double *U, FILE *phifile)		/* DIFFERENT FOR withND */						// function to print the values of the potential and the density in the file tagged as phifile at the given timestep
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

double Int_E(double *U, int i) 		/* DIFFERENT FOR withND */ 						      // Function to calculate the integral of E_h w.r.t. x over the interval I_i = [x_(i-1/2), x_(i+1/2))
{
	int m, j, k;
	double tmp=0., result;
	double ND, eps;																			// declare ND (the value of the doping profile at the given x)
	ND = DopingProfile(i);																	// set ND to the value of the doping profile at i
	eps = EpsilonValue(i);																	// set eps to the value of epsilon at i

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
	result = - tmp*dx*dx*scalev + ND*Gridx((double)i)*dx;
	if(i > a_i)																							// if x > a then there is an extra term to add
	{
		double a_val = (a_i+1)*dx;																		// declare a_val and set it to the value of x at the edge of the ai-th space cell
		result += (NH-NL)*a_val*dx;																		// add (NH-NL)a*dx to result
	}
	if(i > b_i)																							// if x > b then there is an extra term to add
	{
		double b_val = (b_i+1)*dx;																		// declare b_val and set it to the value of x at the edge of the bi-th space cell
		result += (NL-NH)*b_val*dx;																		// add (NL-NH)b*dx to result
	}

	#ifdef Electrons
	result = result/eps - ce*dx;
	#endif	// Electrons
	#ifdef Ions
	result = -result/eps - ce*dx;
	#endif	// Ions
	return result;
}

double Int_E1st(double *U, int i) 	/* DIFFERENT FOR withND */					// \int_i E*(x-x_i)/delta_x dx
{
	int j, k;
	double tmp=0., result;
	double ND, eps;																			// declare ND (the value of the doping profile at the given x)
	ND = DopingProfile(i);																	// set ND to the value of the doping profile at ix
	eps = EpsilonValue(i);																	// set eps to the value of epsilon at ix
	//#pragma omp parallel for reduction(+:tmp)
	for(j=0;j<size_v;j++){
		k=i*size_v + j;
		tmp += U[k*6+0] + U[k*6+5]/4.;
	}
	tmp = tmp*scalev;
	
	#ifdef Electrons
	result = (ND-tmp)*dx*dx/(12.*eps);
	#endif	// Electrons
	#ifdef Ions
	result = (tmp-ND)*dx*dx/(12.*eps);
	#endif	// Ions

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

double Int_E2nd(double *U, int i) 	/* DIFFERENT FOR withND */							// \int_i E* [(x-x_i)/delta_x]^2 dx
{
    int m, j, j1, j2, j3, k;
    double c1=0., c2=0., result;
	double ND, eps;																			// declare ND (the value of the doping profile at the given x)
	ND = DopingProfile(i);																	// set ND to the value of the doping profile at ix
	eps = EpsilonValue(i);																	// set eps to the value of epsilon at i
  
    //cp = computeC_rho(U,i); ce = computePhi_x_0(U);

    for(j=0;j<size_v;j++){
	    k=i*size_v + j;
	    c1 += U[k*6+0] + U[k*6+5]/4.;
	    c2 += U[k*6+1];
    }
    c2 *= dx/2.;				
    
    result = (-cp[i] +scalev*(c1*Gridx(i-0.5) + 0.25*c2))*dx/12. + (ND-scalev*c1)*dx*Gridx((double)i)/12. - scalev*c2*dx/80.;// - ce*dx/12. //BUG: missed -cp

    if(i > a_i)																							// if x > a then there is an extra term to add
	{
		double a_val = (a_i+1)*dx;																		// declare a_val and set it to the value of x at the edge of the ai-th space cell
		result += (NH-NL)*a_val*dx/12.;																	// add (NH-NL)a*dx/12 to result
	}
	if(i > b_i)																							// if x > b then there is an extra term to add
	{
		double b_val = (b_i+1)*dx;																		// declare b_val and set it to the value of x at the edge of the bi-th space cell
		result += (NL-NH)*b_val*dx/12.;																	// add (NL-NH)b*dx/12 to result
	}

	#ifdef Electrons
	result = result/eps - ce*dx/12.;
	#endif	// Electrons
	#ifdef Ions
	result = -result/eps - ce*dx/12.;
	#endif	// Ions

    return result;
}

#endif	/* Doping */
