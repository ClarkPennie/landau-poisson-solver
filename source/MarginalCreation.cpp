/* This is the source file which contains the subroutines necessary for producing marginals of the solutions, suitable for plotting.
 *
 * Functions included: f_marg, PrintMarginalLoc, PrintMarginal
 *
 *  Created on: Nov 15, 2017
 */

double f_marg(double *U, int i, int j1, double x, double v1)
{
	int j2, j3, k0, j2N, k;																		// declare j2, j3 (the indices of the velocity cell in the v2 & v3 directions), k0 (to store the value of i*Nv^3 + j1_Nv^2), j2N (to store the value of j2*N) & k (the index of the cell in U)
	double x_dif, v1_dif, retn;																	// declare x_dif (to store x - x_i), v1_dif (to store v1 - v_j1) & retn (the value of the marginal evaluated at the given x & v1 to be returned at the end

	x_dif = x -  Gridx((double)i);																// set x_dif to x - x_i
	v1_dif = v1 - Gridv((double)j1);															// set v1_dif to v1 - v_j1
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
					fprintf(margfile, "%11.8g  ", x_val);										// in the file tagged as fmarg, print the x coordinate
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
				v1_0 = Gridv((double)j1 - 0.5);													// set v1_0 to the value of v1 at the left edge of the j1-th velocity cell in the v1 direction
				for (nv=0; nv<np; nv++)
				{
					v1_val = v1_0 + nv*ddv;														// set v1_val to v1_0 plus nv increments of width ddv
					fprintf(margfile, "%11.8g  ", v1_val);										// in the file tagged as fmarg, print the v1 coordinate
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
