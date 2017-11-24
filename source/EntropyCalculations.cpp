/* This is the source file which contains the subroutines necessary for calculating the entropy of the
 * solution.
 *
 * Functions included: computeEntropy, computeEntropy_wAvg, computeRelEntropy
 *
 *  Created on: Nov 15, 2017
 */

#include "EntropyCalculations.h"																									// EntropyCalculations.h is where the prototypes for the functions contained in this file are declared

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
