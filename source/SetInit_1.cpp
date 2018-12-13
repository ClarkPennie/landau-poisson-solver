/* This is the source file which contains the subroutines necessary for setting the initial conditions, as well as
 * projecting to the function required as the initial condition to the spectral method for the collision problem
 * resulting from time-splitting, as well as the function for setting the trapezoidal weights.
 *
 * Functions included: trapezoidalRule, f_TS, f_2Gauss, Mw, Mw_x, f_2H, SetInit_LD, SetInit_4H, SetInit_2H,
 * setInit_spectral
 *
 */

//static double vt[4] = {-0.3399810435848562648026658,0.3399810435848562648026658,-0.8611363115940525752239465,0.8611363115940525752239465};
//static double wt[4] = {0.6521451548625461426269361,0.6521451548625461426269361,0.3478548451374538573730639,0.3478548451374538573730639};
//double wt[5]={0.5688888888888889, 0.4786286704993665, 0.4786286704993665,0.2369268850561891, 0.2369268850561891};				// weights for Gaussian quadrature
//double vt[5]={0., -0.5384693101056831,0.5384693101056831,-0.9061798459386640,0.9061798459386640};								// node values for Gaussian quadrature over the interval [-1,1]
	
#include "SetInit_1.h"																							// SetInit_1.h is where the prototypes for the functions contained in this file are declared

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

double Mw(double v1, double v2, double v3)																		// function to return a Maxwellian in v-space
{
	double r2, T, retn;																							// declare r2 (the squared magnitude of v for the given v = (v1,v2,v3)), T (the temperature) and retn (the value of the Maxwellian evaluated at (v1,v2,v3) to be returned)
	T=0.4; 																										// set T to 0.9 (when testing the nonlinear damping, T was chosen too small that the "effective grid" is  not fine enough)
	r2=v1*v1+v2*v2+v3*v3;																						// calculate r2 for the given (v1,v2,v3)
	retn = exp(-r2/(2*T))/(2*PI*T*sqrt(2*T*PI));																// calculate the value of the Maxwellian at the given point (v1,v2,v3)
	return retn;																								// return the result
}

double Mw_x(double x)																							// function to return a Maxwellian in x-space
{
	double r2, T, retn;																							// declare r2 (the squared magnitude of the given x), T (the temperature) and retn (the value of the Maxwellian evaluated at (v1,v2,v3) to be returned)
	T=0.4; 																										// set T to 0.9 (when testing the nonlinear damping, T was chosen too small that the "effective grid" is  not fine enough)
	r2=x*x;																										// calculate r2 for the given x
	retn = exp(-r2/(2*T))/(sqrt(2*T*PI));																		// calculate the value of the Maxwellian at the given point x
	return retn;																								// return the result
}

double f_2H(double x)																							// function to calculate the value of f_DH(x) = sin^2((x-Lx/2)^2/Lx)*exp(-((x-Lx/2)^2)/2T)/(c*sqrt(2*pi*T)) (to create a double hump function, where c = 0.0145254 normalises mass to 1)
{
	double r2, MW, retn;																						// declare r2 (the squared magnitude of the given x), MW (the value of the Maxwellian evaluated at x) and retn (the value of f_DH(x))
	r2 = (x-Lx/2)*(x-Lx/2);																						// calculate r2 for the given x
	MW = Mw_x(x-Lx/2);																							// calculate the Maxwellian in x-space evaluated at the given x
	retn = sin(r2/Lx)*sin(r2/Lx)*MW/0.0145254;																	// calculate the value of f_DH(x)
	return retn;																								// return the result
}

void SetInit_LD(double *U)																						// function to calculate the DG coefficients for the initial condition for Landau Damping
{
    int i, j1, j2, j3, k, m1,m2,m3,nt=5;																		// declare i (to represent cell i in x-space), j1, j2, j3 (to represent cell (j1,j2,j3) in v-space), k (the index of cell (i,j1,j2,j3) in U), m1, m2, m3 (counters for the Gaussian quadrature in 3D) & nt (the number of points in the quadrature)
    double a=A_amp, c=k_wave;																					// declare a (the amplitude of cosine wave) and set it to A_amp & c (the frequency of the cosine wave) and set it to k_wave
    double tp, tp0, tp5, tmp0, tmp1, tmp2, tmp3, tmp4;															// declare tp, tp0, tmp0, tmp1, tmp2, tmp3, tmp4 (temporary values while calculating the quadrature for the integral w.r.t. v)
    //#pragma omp parallel for private(k,j1,j2,j3,i,tmp0, tmp1, tmp2, tmp3, tmp4, tp0, tp5, tp) shared(U)
    for(j1=0;j1<Nv;j1++)																						// loop through all the velocity cells
    {
    	for(j2=0;j2<Nv;j2++)
    	{
    		for(j3=0;j3<Nv;j3++)
    		{
    			tmp0=0.; tmp1=0.; tmp2=0.; tmp3=0.; tmp4=0.;													// initialise tmp0, tmp1, tmp2, tmp3 & tmp4 at 0 for a new quadrature integral to calculate int_Kj Mw(v)*phi_(6k+l) dv, for l = 0, 2, 4, 5, 6
    			for(m1=0;m1<nt;m1++)																			// loop through the quadrature sum
    			{
    				for(m2=0;m2<nt;m2++)
    				{
    					for(m3=0;m3<nt;m3++)
    					{
							if(Damping)
							{
								tp = wt[m1]*wt[m2]*wt[m3]*Mw(Gridv((double)j1)+0.5*dv*vt[m1], Gridv((double)j2)+0.5*dv*vt[m2], Gridv((double)j3)+0.5*dv*vt[m3]);		// calculate w_m1*w_m2*w_m3*Mw(v_m1,v_m2,v_m3), which appears in all quadrature integral approximations
							}

							if(TwoStream)
							{
								tp = wt[m1]*wt[m2]*wt[m3]*f_2Gauss(Gridv((double)j1)+0.5*dv*vt[m1], Gridv((double)j2)+0.5*dv*vt[m2], Gridv((double)j3)+0.5*dv*vt[m3]);
							}

							tmp0 += tp;																																// add tp to tmp0 (for the integral int_Kj Mw(v)*phi_6k dv)
							tmp1 += tp*0.5*vt[m1];																													// add tp*v_m1/2 to tmp1 (for the integral int_Kj Mw(v)*phi_(6k+2) dv)
							tmp2 += tp*0.5*vt[m2];																													// add tp*v_m2/2 to tmp2 (for the integral int_Kj Mw(v)*phi_(6k+3) dv)
							tmp3 += tp*0.5*vt[m3];																													// add tp*v_m3/2 to tmp3 (for the integral int_Kj Mw(v)*phi_(6k+4) dv)
							tmp4 += tp*0.25*(vt[m1]*vt[m1] + vt[m2]*vt[m2]+ vt[m3]*vt[m3]);																			// add tp*((v_m1/2)^2 + (v_m2/2)^2 + (v_m3/2)^2) to tmp4 (for the integral int_Kj Mw(v)*phi_(6k+5) dv)
    					}
    				}
    			}
    			tmp0 = tmp0*0.5*0.5*0.5; tmp1 = tmp1*0.5*0.5*0.5; tmp2 = tmp2*0.5*0.5*0.5; tmp3 = tmp3*0.5*0.5*0.5; tmp4 = tmp4*0.5*0.5*0.5;						// multiply tmp0, tmp1, tmp2, tmp3 & tmp4 by (1/2)^3 to represent the fact that quadrature isn't done over [-1, 1] (should also multiply by dv^3 but this cancels with 1/dv^3 later)
    			for(i=0;i<Nx;i++)																																	// loop through the space cells
    			{
    				k=i*size_v + (j1*Nv*Nv + j2*Nv + j3);																											// calculate the index of cell (i,j1,j2,j3) in U
    				tp0 = (dx + (sin(c*Gridx(i+0.5)) - sin(c*Gridx(i-0.5)))*a/c)*tmp0/dx;																			// calculate b_6k = (int_Ii (1 + Acos(kx)) dx)*(int_Kj Mw(v)*phi_6k(v) dv) (NOTE: int_(Omega_i) (1 + a*cos(c*x)) dx = dx + (sin(c*x_(i+0.5)) - sin(c*x_(i+0.5)))*(a/c))
    				tp5 = (dx + (sin(c*Gridx(i+0.5)) - sin(c*Gridx(i-0.5)))*a/c)*tmp4/dx;																			// calculate b_(6k+5) = (int_Ii (1 + Acos(kx)) dx)*(int_Kj Mw(v)*phi_(6k+5)(v) dv) (NOTE: int_(Omega_i) (1 + a*cos(c*x)) dx = dx + (sin(c*x_(i+0.5)) - sin(c*x_(i+0.5)))*(a/c))
    				U[k*6+0] = 19*tp0/4. - 15*tp5;																													// calculate the coefficient U[6k]
    				U[k*6+5] = 60*tp5 - 15*tp0;																														// calculate the coefficient U[6k+5]

    				U[k*6+1] = (0.5*(sin(c*Gridx(i+0.5)) + sin(c*Gridx(i-0.5))) + (cos(c*Gridx(i+0.5)) - cos(c*Gridx(i-0.5)))/(c*dx))*(a/c)*tmp0*12./dx; 			// calculate the coefficient U[6k+1] = 12*b_(6k+1) = 12*(int_Ii (1 + Acos(kx))*phi_(6k+1)(x) dx)*(int_Kj Mw(v) dv) (NOTE: int_(Omega_i) (1 + a*cos(c*x))*phi_(6k+1)(x) dx = (0.5*(sin(c*x_(i+0.5)) + sin(c*x_(i+0.5))) + (cos(c*x_(i+0.5)) - cos(c*x_(i+0.5)))/(c*dx))*(a/c))
    				U[k*6+2] = (dx + (sin(c*Gridx(i+0.5)) - sin(c*Gridx(i-0.5)))*a/c)*tmp1*12/dx;																	// calculate the coefficient U[6k+2] = 12*b_(6k+2) = 12*(int_Ii (1 + Acos(kx)) dx)*(int_Kj Mw(v)*phi_(6k+2)(v) dv) (NOTE: int_(Omega_i) (1 + a*cos(c*x)) dx = dx + (sin(c*x_(i+0.5)) - sin(c*x_(i+0.5)))*(a/c))
    				U[k*6+3] = (dx + (sin(c*Gridx(i+0.5)) - sin(c*Gridx(i-0.5)))*a/c)*tmp2*12/dx;																	// calculate the coefficient U[6k+3] = 12*b_(6k+3) = 12*(int_Ii (1 + Acos(kx)) dx)*(int_Kj Mw(v)*phi_(6k+3)(v) dv) (NOTE: int_(Omega_i) (1 + a*cos(c*x)) dx = dx + (sin(c*x_(i+0.5)) - sin(c*x_(i+0.5)))*(a/c))
    				U[k*6+4] = (dx + (sin(c*Gridx(i+0.5)) - sin(c*Gridx(i-0.5)))*a/c)*tmp3*12/dx;																	// calculate the coefficient U[6k+4] = 12*b_(6k+4) = 12*(int_Ii (1 + Acos(kx)) dx)*(int_Kj Mw(v)*phi_(6k+4)(v) dv) (NOTE: int_(Omega_i) (1 + a*cos(c*x)) dx = dx + (sin(c*x_(i+0.5)) - sin(c*x_(i+0.5)))*(a/c))
    			}
    		}
    	}
	}
}

void SetInit_4H(double *U)																							// function to calculate the DG coefficients for the initial condition with four humps, found by adding four Maxwellians
{
    int i, j1, j2, j3, k, m,m1,m2,m3,nt=5, p;																		// declare i (to represent cell i in x-space), j1, j2, j3 (to represent cell (j1,j2,j3) in v-space), k (the index of cell (i,j1,j2,j3) in U), m (counter for the Gaussian quadrature in x-space), m1, m2, m3 (counters for the Gaussian quadrature in v-space), nt (the number of points in the quadrature) & p (to loop through the four Maxwellians)
    double tp, tpx, tp0, tp5, tmpx0, tmpx1, tmp0, tmp1, tmp2, tmp3, tmp4;											// declare tp, tpx, tp0, tmpx0, tmpx1, tmp0, tmp1, tmp2, tmp3, tmp4 (temporary values while calculating the quadrature for the integral w.r.t. v)
    //#pragma omp parallel for private(k,j1,j2,j3,i,tmp0, tmp1, tmp2, tmp3, tmp4, tp0, tp5, tp) shared(U)
    for(p=0;p<4;p++)																								// loop through the four Maxwellians
    {
    	for(j1=0;j1<Nv;j1++)																						// loop through all the velocity cells
    	{
    		for(j2=0;j2<Nv;j2++)
    		{
    			for(j3=0;j3<Nv;j3++)
    			{
    				tmp0=0.; tmp1=0.; tmp2=0.; tmp3=0.; tmp4=0.;													// initialise tmp0, tmp1, tmp2, tmp3 & tmp4 at 0 for a new quadrature integral to calculate int_Kj Mw(v)*phi_(6k+l)(v) dv, for l = 0, 2, 4, 5, 6
    				for(m1=0;m1<nt;m1++)																			// loop through the quadrature sum
    				{
    					for(m2=0;m2<nt;m2++)
    					{
    						for(m3=0;m3<nt;m3++)
    						{
    							tp = wt[m1]*wt[m2]*wt[m3]*Mw(Gridv((double)j1)+0.5*dv*vt[m1] + pow(-1,p), Gridv((double)j2)+0.5*dv*vt[m2] + pow(-1,p), Gridv((double)j3)+0.5*dv*vt[m3] + pow(-1,p));		// calculate w_m1*w_m2*w_m3*Mw(v_m1+(-1)^p,v_m2+(-1)^p,v_m3+(-1)^p), a Maxwellian shifted to center at v_j = (-1)^p, which appears in all quadrature integral approximations

    							tmp0 += tp;																			// add tp to tmp0 (for the integral int_Kj Mw(v) dv)
    							tmp1 += tp*0.5*vt[m1];																// add tp*v_m1/2 to tmp1 (for the integral int_Kj Mw(v)*phi_(6k+2)(v) dv)
    							tmp2 += tp*0.5*vt[m2];																// add tp*v_m2/2 to tmp2 (for the integral int_Kj Mw(v)*phi_(6k+3)(v) dv)
    							tmp3 += tp*0.5*vt[m3];																// add tp*v_m3/2 to tmp3 (for the integral int_Kj Mw(v)*phi_(6k+4)(v) dv)
    							tmp4 += tp*0.25*(vt[m1]*vt[m1] + vt[m2]*vt[m2]+ vt[m3]*vt[m3]);						// add tp*((v_m1/2)^2 + (v_m2/2)^2 + (v_m3/2)^2) to tmp4 (for the integral int_Kj Mw(v)*phi_(6k+5)(v) dv)
    						}
    					}
    				}
    				tmp0 = tmp0*0.5*0.5*0.5; tmp1 = tmp1*0.5*0.5*0.5; tmp2 = tmp2*0.5*0.5*0.5; tmp3 = tmp3*0.5*0.5*0.5; tmp4 = tmp4*0.5*0.5*0.5;						// multiply tmp0, tmp1, tmp2, tmp3 & tmp4 by (1/2)^3 to represent the fact that quadrature isn't done over [-1, 1] (should also multiply by dv^3 but this cancels with 1/dv^3 later)
    				for(i=0;i<Nx;i++)																				// loop through the space cells
    				{
    					tmpx0 = 0.; tmpx1 = 0.;																		// initialise tmpx0 & tmpx1 at 0 for a new quadrature integral to calculate int_Ii f_DH(x)*phi_(6k+l)(x) dx, for l = 0, 1
    					for(m = 0; m < nt; m++)																		// loop through the quadrature sum
    					{
    						tpx = wt[m]*Mw_x(Gridx((double) i)+0.5*dx*vt[m] - Lx/2 + pow(-1,((int)(p/2))));	// calculate w_m*Mw_x(x-Lx/2+(-1)^floor(p/2)), a Maxwellian shifted to center at x = Lx/2 - (-1)^floor(p/2), which appears in both quadrature integral approximations
    						tmpx0 += tpx;																			// add tpx to tmpx0 (for the integral int_Ii f_DH(x) dx)
    						tmpx1 += tpx*0.5*vt[m];																	// add tpx*x_m/2 to tmpx1 (for the integral int_Ii f_DH(x)*phi_(6k+1)(x) dx)
    					}
    					tmpx0 = tmpx0*0.5; tmpx1 = tmpx1*0.5;														// multiply tmpx0 & tmpx1 by 1/2 to represent the fact that quadrature isn't done over [-1, 1] (should also multiply by dx but this cancels with 1/dx later)
    					k=i*size_v + (j1*Nv*Nv + j2*Nv + j3);														// calculate the index of cell (i,j1,j2,j3) in U

    					tp0 = tmpx0*tmp0;																			// calculate b_6k = (int_Ii f_DH(x) dx)*(int_Kj Mw(v) dv)
    					tp5 = tmpx0*tmp4;																			// calculate b_(6k+5) = (int_Ii f_DH(x) dx)*(int_Kj Mw(v)*phi_(6k+5)(v) dv)

    					if(p==0)
    					{
    						U[k*6+0] = 19*tp0/4. - 15*tp5;															// calculate the coefficient U[6k]
    						U[k*6+5] = 60*tp5 - 15*tp0;																// calculate the coefficient U[6k+5]

    						U[k*6+1] = tmpx1*tmp0*12; 																// calculate the coefficient U[6k+1] = 12*b_(6k+1) = 12*(int_Ii f_DH(x)*phi_(6k+1)(x) dx)*(int_Kj Mw(v) dv)
    						U[k*6+2] = tmpx0*tmp1*12;																// calculate the coefficient U[6k+2] = 12*b_(6k+2) = 12*(int_Ii f_DH(x) dx)*(int_Kj Mw(v)*phi_(6k+2)(v) dv)
    						U[k*6+3] = tmpx0*tmp2*12;																// calculate the coefficient U[6k+3] = 12*b_(6k+3) = 12*(int_Ii f_DH(x) dx)*(int_Kj Mw(v)*phi_(6k+3)(v) dv)
    						U[k*6+4] = tmpx0*tmp3*12;																// calculate the coefficient U[6k+4] = 12*b_(6k+4) = 12*(int_Ii f_DH(x) dx)*(int_Kj Mw(v)*phi_(6k+4)(v) dv)
    					}
    					else
    					{
    						U[k*6+0] += 19*tp0/4. - 15*tp5;															// calculate the coefficient U[6k]
    						U[k*6+5] += 60*tp5 - 15*tp0;															// calculate the coefficient U[6k+5]

    						U[k*6+1] += tmpx1*tmp0*12; 																// calculate the coefficient U[6k+1] = 12*b_(6k+1) = 12*(int_Ii f_DH(x)*phi_(6k+1)(x) dx)*(int_Kj Mw(v) dv)
    						U[k*6+2] += tmpx0*tmp1*12;																// calculate the coefficient U[6k+2] = 12*b_(6k+2) = 12*(int_Ii f_DH(x) dx)*(int_Kj Mw(v)*phi_(6k+2)(v) dv)
    						U[k*6+3] += tmpx0*tmp2*12;																// calculate the coefficient U[6k+3] = 12*b_(6k+3) = 12*(int_Ii f_DH(x) dx)*(int_Kj Mw(v)*phi_(6k+3)(v) dv)
    						U[k*6+4] += tmpx0*tmp3*12;																// calculate the coefficient U[6k+4] = 12*b_(6k+4) = 12*(int_Ii f_DH(x) dx)*(int_Kj Mw(v)*phi_(6k+4)(v) dv)
    					}
    				}
    			}
    		}
    	}
    }
    for(k=0;k<6*size;k++)																							// loop through all entries of U
    {
    	U[k] = U[k]/4;																								// divide the kth entry of U by 4 since 4 Maxwellians were added together
    }
}

void SetInit_2H(double *U)																							// function to calculate the DG coefficients for the initial condition with two humps
{
    int i, j1, j2, j3, k, m,m1,m2,m3,nt=5;																			// declare i (to represent cell i in x-space), j1, j2, j3 (to represent cell (j1,j2,j3) in v-space), k (the index of cell (i,j1,j2,j3) in U), m (counter for the Gaussian quadrature in x-space), m1, m2, m3 (counters for the Gaussian quadrature in v-space) & nt (the number of points in the quadrature)
    double tp, tpx, tp0, tp5, tmpx0, tmpx1, tmp0, tmp1, tmp2, tmp3, tmp4;											// declare tp, tpx, tp0, tmpx0, tmpx1, tmp0, tmp1, tmp2, tmp3, tmp4 (temporary values while calculating the quadrature for the integral w.r.t. v)
    //#pragma omp parallel for private(k,j1,j2,j3,i,tmp0, tmp1, tmp2, tmp3, tmp4, tp0, tp5, tp) shared(U)
    for(j1=0;j1<Nv;j1++)																							// loop through all the velocity cells
    {
    	for(j2=0;j2<Nv;j2++)
    	{
    		for(j3=0;j3<Nv;j3++)
    		{
    			tmp0=0.; tmp1=0.; tmp2=0.; tmp3=0.; tmp4=0.;														// initialise tmp0, tmp1, tmp2, tmp3 & tmp4 at 0 for a new quadrature integral to calculate int_Kj Mw(v)*phi_(6k+l)(v) dv, for l = 0, 2, 4, 5, 6
    			for(m1=0;m1<nt;m1++)																				// loop through the quadrature sum
    			{
    				for(m2=0;m2<nt;m2++)
    				{
    					for(m3=0;m3<nt;m3++)
    					{
							tp = wt[m1]*wt[m2]*wt[m3]*Mw(Gridv((double)j1)+0.5*dv*vt[m1], Gridv((double)j2)+0.5*dv*vt[m2], Gridv((double)j3)+0.5*dv*vt[m3]);		// calculate w_m1*w_m2*w_m3*Mw(v_m1,v_m2,v_m3), which appears in all quadrature integral approximations

							tmp0 += tp;																				// add tp to tmp0 (for the integral int_Kj Mw(v) dv)
							tmp1 += tp*0.5*vt[m1];																	// add tp*v_m1/2 to tmp1 (for the integral int_Kj Mw(v)*phi_(6k+2)(v) dv)
							tmp2 += tp*0.5*vt[m2];																	// add tp*v_m2/2 to tmp2 (for the integral int_Kj Mw(v)*phi_(6k+3)(v) dv)
							tmp3 += tp*0.5*vt[m3];																	// add tp*v_m3/2 to tmp3 (for the integral int_Kj Mw(v)*phi_(6k+4)(v) dv)
							tmp4 += tp*0.25*(vt[m1]*vt[m1] + vt[m2]*vt[m2]+ vt[m3]*vt[m3]);							// add tp*((v_m1/2)^2 + (v_m2/2)^2 + (v_m3/2)^2) to tmp4 (for the integral int_Kj Mw(v)*phi_(6k+5)(v) dv)
    					}
    				}
    			}
    			tmp0 = tmp0*0.5*0.5*0.5; tmp1 = tmp1*0.5*0.5*0.5; tmp2 = tmp2*0.5*0.5*0.5; tmp3 = tmp3*0.5*0.5*0.5; tmp4 = tmp4*0.5*0.5*0.5;						// multiply tmp0, tmp1, tmp2, tmp3 & tmp4 by (1/2)^3 to represent the fact that quadrature isn't done over [-1, 1] (should also multiply by dv^3 but this cancels with 1/dv^3 later)
    			for(i=0;i<Nx;i++)																					// loop through the space cells
    			{
    				tmpx0 = 0.; tmpx1 = 0.;																			// initialise tmpx0 & tmpx1 at 0 for a new quadrature integral to calculate int_Ii f_DH(x)*phi_(6k+l)(x) dx, for l = 0, 1
    				for(m = 0; m < nt; m++)																			// loop through the quadrature sum
    				{
    					tpx = wt[m]*f_2H(Gridx((double) i)+0.5*dx*vt[m]);											// calculate w_m*f_DH(x), which appears in both quadrature integral approximations
    					tmpx0 += tpx;																				// add tpx to tmpx0 (for the integral int_Ii f_DH(x) dx)
    					tmpx1 += tpx*0.5*vt[m];																		// add tpx*x_m/2 to tmpx1 (for the integral int_Ii f_DH(x)*phi_(6k+1)(x) dx)
    				}
    				tmpx0 = tmpx0*0.5; tmpx1 = tmpx1*0.5;															// multiply tmpx0 & tmpx1 by 1/2 to represent the fact that quadrature isn't done over [-1, 1] (should also multiply by dx but this cancels with 1/dx later)
    				k=i*size_v + (j1*Nv*Nv + j2*Nv + j3);															// calculate the index of cell (i,j1,j2,j3) in U

    				tp0 = tmpx0*tmp0;																				// calculate b_6k = (int_Ii f_DH(x) dx)*(int_Kj Mw(v) dv)
    				tp5 = tmpx0*tmp4;																				// calculate b_(6k+5) = (int_Ii f_DH(x) dx)*(int_Kj Mw(v)*phi_(6k+5)(v) dv)
    				U[k*6+0] = 19*tp0/4. - 15*tp5;																	// calculate the coefficient U[6k]
    				U[k*6+5] = 60*tp5 - 15*tp0;																		// calculate the coefficient U[6k+5]

    				U[k*6+1] = tmpx1*tmp0*12; 																		// calculate the coefficient U[6k+1] = 12*b_(6k+1) = 12*(int_Ii f_DH(x)*phi_(6k+1)(x) dx)*(int_Kj Mw(v) dv)
    				U[k*6+2] = tmpx0*tmp1*12;																		// calculate the coefficient U[6k+2] = 12*b_(6k+2) = 12*(int_Ii f_DH(x) dx)*(int_Kj Mw(v)*phi_(6k+2)(v) dv)
    				U[k*6+3] = tmpx0*tmp2*12;																		// calculate the coefficient U[6k+3] = 12*b_(6k+3) = 12*(int_Ii f_DH(x) dx)*(int_Kj Mw(v)*phi_(6k+3)(v) dv)
    				U[k*6+4] = tmpx0*tmp3*12;																		// calculate the coefficient U[6k+4] = 12*b_(6k+4) = 12*(int_Ii f_DH(x) dx)*(int_Kj Mw(v)*phi_(6k+4)(v) dv)
    			}
    		}
    	}
	}
}

#ifdef UseMPI
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
