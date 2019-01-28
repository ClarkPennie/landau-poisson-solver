/* This is the source file which contains the subroutines necessary for checking where the solution
 * loses positivity.
 *
 * Functions included: computeCellAvg, FindNegVals, FindNegVals_old, CheckNegVals
 *
 *  Created on: Nov 15, 2017
 */

#include "NegativityChecks.h"																										// NegativityChecks.h is where the prototypes for the functions contained in this file are declared

/*
double computeCellAvg(double *U, int i, int j1, int j2, int j3)
{
	if(Homogeneous)
	{
		return computeCellAvg_Homo(U, j1, j2, j3);
	}
	else
	{
		return computeCellAvg_Inhomo(U, i, j1, j2, j3);
	}
}
*/

double computeCellAvg_Inhomo(double *U, int i, int j1, int j2, int j3)																		// function to calculate the average value of the approximate function f (with DG coefficients in U) on the cell I_i x K_(j1,j2,j3), namely (1/cell_volume)*int_(I_i x K_(j1,j2,j3)) f dxdv = (1/(dx*dv^3))*int_(I_i x K_(j1,j2,j3)) f dxdv
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

double computeCellAvg_Homo(double *U, int j1, int j2, int j3)																			// function to calculate the average value of the approximate function f (with DG coefficients in U) on the cell K_(j1,j2,j3), namely (1/cell_volume)*int_(K_(j1,j2,j3)) f dv = (1/(dv^3))*int_(K_(j1,j2,j3)) f dv
{
	int k, nv1, nv2, nv3;																											// declare k (the location of the given cell in U) & nv1, nv2, nv3 (counters for the velocity cell)
	double v1_0, v2_0, v3_0, v1_val, v2_val, v3_val, f_val, avg;																	// declare v1_0, v2_0, v3_0 (to store the values in the middle of the current cell), v1_val, v2_val, v3_val (to store the v values to be evaluated at), f_val (to store the value of the function evaluated at the current v values) & avg (to store the average to be returned)
	k = Nv*(Nv*j1 + j2) + j3;																										// set k to j1*Nv^2 + j2*Nv + j3
	v1_0 = Gridv((double)j1);																										// set v1_0 to the value of v1 at the center of the j1th velocity cell in that direction
	v2_0 = Gridv((double)j2);																										// set v2_0 to the value of v2 at the center of the j2th velocity cell in that direction
	v3_0 = Gridv((double)j3);																										// set v3_0 to the value of v3 at the center of the j3th velocity cell in that direction
	avg = 0;																														// set avg to 0 originally

	//#pragma omp parallel for private(nx,nv1,nv2,nv3,x_val,v1_val,v2_val,v3_val,f_val) shared(U,dv,dx,vt,wt,x_0,v1_0,v2_0,v3_0,k) reduction(+:avg)
	for(nv1=0;nv1<5;nv1++)																											// loop through the five quadrature points in the v1 direction of the cell
	{
		v1_val = v1_0 + 0.5*vt[nv1]*dv;																								// set x_val to the nv1-th quadrature point in the cell
		for(nv2=0;nv2<5;nv2++)																										// loop through the five quadrature points in the v2 direction of the cell
		{
			v2_val = v2_0 + 0.5*vt[nv2]*dv;																							// set x_val to the nv2-th quadrature point in the cell
			for(nv3=0;nv3<5;nv3++)																									// loop through the five quadrature points in the v3 direction of the cell
			{
				v3_val = v3_0 + 0.5*vt[nv3]*dv;																						// set x_val to the nv3-th quadrature point in the cell
				f_val = U[k*6+0] + U[k*6+2]*(v1_val-v1_0)/dv +
						U[k*6+3]*(v2_val-v2_0)/dv + U[k*6+4]*(v3_val-v3_0)/dv +
						U[k*6+5]*(((v1_val-v1_0)/dv)*((v1_val-v1_0)/dv)+((v2_val-v2_0)/dv)*((v2_val-v2_0)/dv)
								+((v3_val-v3_0)/dv)*((v3_val-v3_0)/dv));															// set f_val to the evaluation of the approximation at (x_val,v1_val,v2_val,v3_val)
				avg += wt[nv1]*wt[nv2]*wt[nv3]*f_val;																				// add f(x_val,v1_val,v2_val,v3_val) times the quadrature weights corresponding to the current x & v values to the current quadrature value
			}
		}
	}
	avg = avg*0.5*0.5*0.5;																											// multiply the quadrature result by 0.5*0.5*0.5*0.5, since each quadrature is done over intervals of width dv (three times) or dx (once) instead of the standard interval [-1,1] (should be multiplied by dx*dv^3 also but then, to find the average, need to divide by dx*dv^3)
	return avg;																														// return the average value of f on the cell
}

void FindNegVals(double *U, int *NegVals, double *AvgVals)
{
	if(Homogeneous)
	{
		FindNegVals_Homo(U, NegVals, AvgVals);
	}
	else
	{
		FindNegVals_Inhomo(U, NegVals, AvgVals);
	}
}

void FindNegVals_Inhomo(double *U, int *NegVals, double *AvgVals)																							// function to find out the cells in which the approximation from U is negative on average and stores the cell locations in NegVals
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
					f_avg = computeCellAvg_Inhomo(U,i,j1,j2,j3);																			// calculate the average value of the approximate solution in the current cell and set it to f_avg
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

void FindNegVals_Homo(double *U, int *NegVals, double *AvgVals)																							// function to find out the cells in which the approximation from U is negative on average and stores the cell locations in NegVals
{
	int j1, j2, j3, j1NN, j2N, k, nv1, nv2, nv3;																					// declare j1, j2, j3 (the indices of the velocity cell), j1NN (to store j1*Nv^2), j2N (to store j2*Nv), k (the location of the given cell in U) & nv1, nv2, nv3 (counters for the velocity cell)
	double v1_val, v2_val, v3_val, f_avg;																							// declare v1_val, v2_val, v3_val (to store the v values to be evaluated at) & f_val (to store the value of the function evaluated at the current v values)
	#pragma omp parallel for private(j1,j2,j3,k,j1NN,j2N,f_avg) shared(U,NegVals,AvgVals,size_v,Nv)
	for(j1=0;j1<Nv;j1++)																											// loop through the velocity cells in the v1 direction
	{
		j1NN = j1*Nv*Nv;																											// set j1NN to j1*Nv^2
		for(j2=0;j2<Nv;j2++)																										// loop through the velocity cells in the v2 direction
		{
			j2N = j2*Nv;																											// set j2N to j2*Nv
			for(j3=0;j3<Nv;j3++)																									// loop through the velocity cells in the v3 direction
			{
				k = j1NN + j2N + j3;																								// set k to j1*Nv^2 + j2*Nv + j3
				NegVals[k] = 0;																										// set NegVals[k] to 0, which assumes at first that there will be no negative values in the kth cell
				f_avg = computeCellAvg_Homo(U,j1,j2,j3);																					// calculate the average value of the approximate solution in the current cell and set it to f_avg
				AvgVals[k] = f_avg;																									// store f_avg in AvgVals[k]
				if(f_avg < 0)																										// check if this value was negative
				{
					NegVals[k] = 1;																									// if so, set NegVals[k] to 1 to indicate that there was a negative value in this cell
				}
			}
		}
	}
}

/*
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
*/

void CheckNegVals(double *U, int *NegVals, double *AvgVals)
{
	if(Homogeneous)
	{
		CheckNegVals_Homo(U, NegVals, AvgVals);
	}
	else
	{
		CheckNegVals_Inhomo(U, NegVals, AvgVals);
	}
}

void CheckNegVals_Inhomo(double *U, int *NegVals, double *AvgVals)																							// function to find out the cells in which the approximation from U turns negative and stores the cell locations in NegVals
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

void CheckNegVals_Homo(double *U, int *NegVals, double *AvgVals)																							// function to find out the cells in which the approximation from U turns negative and stores the cell locations in NegVals
{
	int j1, j2, j3, j1NN, j2N, k;																								// declare j1, j2, j3 (the indices of the velocity cell), j1NN (to store j1*Nv^2), j2N (to store j2*Nv), k (the location of the given cell in U) & nv1, nv2, nv3 (counters for the velocity cell)
	double f_avg;																													// declare f_val (to store the value of the function evaluated at the current x & v values)
	for(j1=0;j1<Nv;j1++)																											// loop through the velocity cells in the v1 direction
	{
		j1NN = j1*Nv*Nv;																											// set j1NN to j1*Nv^2
		for(j2=0;j2<Nv;j2++)																										// loop through the velocity cells in the v2 direction
		{
			j2N = j2*Nv;																											// set j2N to j2*Nv
			for(j3=0;j3<Nv;j3++)																									// loop through the velocity cells in the v3 direction
			{
				k = j1NN + j2N + j3;																								// set k to j1*Nv^2 + j2*Nv + j3
				printf("%d ~ (%d, %d, %d) f_avg = %g \n", NegVals[k], j1, j2, j3, AvgVals[k]);
				printf(" \n");
			}
		}
	}
}


