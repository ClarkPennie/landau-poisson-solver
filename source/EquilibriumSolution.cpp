/* This is the source file which contains the subroutines necessary for constructing an equilibrium solution (necessary for the relative entropy calculations).
 *
 * Functions included: ExportRhoQuadVals, ComputeEquiVals, PrintEquiVals
 *
 *  Created on: Nov 15, 2017
 */

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
	/*DEGUG TEST
	for(i=0;i<5*Nx;i++)
	{
		printf("%g ", rho_vals[i]);
	}i
	*/
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
