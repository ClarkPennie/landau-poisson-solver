/* This is the source file which contains the subroutines necessary for calculating the moments of the
 * solution.
 *
 * Functions included: computeMass, computeMomentum, computeKiE, computeKiEratio, computeEleE
 *
 *  Created on: Nov 15, 2017
 */

#include "MomentCalculations.h"																			// MomentCalculations.h is where the prototypes for the functions contained in this file are declared

double computeMass(double *U)
{
	if(Homogeneous)
	{
		return computeMass_Homo(U);
	}
	else
	{
		return computeMass_Inhomo(U);
	}
}

double computeMass_Inhomo(double *U)
{
	if(MeshRefinement)
	{
		return computeMass_Inhomo_Refined(U);
	}
	else
	{
		return computeMass_Inhomo_Uniform(U);
	}
}

double computeMass_Inhomo_Uniform(double *U)
{
  int k;
  double tmp=0.;
  #pragma omp parallel for shared(U) reduction(+:tmp)
  for(k=0;k<Nx*size_v;k++) tmp += U[k*6+0] + U[k*6+5]/4.;

  return tmp*dx*scalev;
}

double computeMass_Inhomo_Refined(double *U)
{
	int i, k, k0, k_v;
	double dx_val;
	double tmp=0.;
	#pragma omp parallel for private(i, k_v, k0, k, dx_val) shared(U) reduction(+:tmp)
	for(i=0; i<Nx; i++)
	{
		dx_val = dx_value(i);
		k0 = i*size_v;

		for(k_v=0; k_v<size_v; k_v++)
		{
			k = k0 + k_v;
			tmp += dx_val*(U[k*6+0] + U[k*6+5]/4.);
		}
	}
	return tmp*scalev;
}

double computeMass_Homo(double *U)
{
  int k;
  double tmp=0.;
  #pragma omp parallel for shared(U) reduction(+:tmp)
  for(k=0;k<size_v;k++) tmp += U[k*6+0] + U[k*6+5]/4.;

  return tmp*scalev;
}

double computeMass_in_x(double *U, int i, double x)
{
	int k, k0, k_v;
	double dx_val, phi1_val;
	double tmp=0.;

	dx_val = dx_value(i);
	phi1_val = (x -  Gridx((double)i))/dx_val;

	k0 = i*Nv*Nv*Nv;

	#pragma omp parallel for private(k_v, k) shared(U) reduction(+:tmp)
	for(k_v=0;k_v<size_v;k_v++)
	{
	  k = k0 + k_v;
	  tmp += U[k*6+0] + U[k*6+1]*phi1_val + U[k*6+5]/4.;
	}

	return tmp*scalev;
}

void computeMomentum(double *U, double *a)
{
	if(Homogeneous)
	{
		return computeMomentum_Homo(U, a);
	}
	else
	{
		return computeMomentum_Inhomo(U, a);
	}
}

void computeMomentum_Inhomo(double *U, double *a)
{
	if(MeshRefinement)
	{
		return computeMomentum_Inhomo_Refined(U, a);
	}
	else
	{
		return computeMomentum_Inhomo_Uniform(U, a);
	}
}

void computeMomentum_Inhomo_Refined(double *U, double *a)
{
	int k, k0, k_v, i, j1,j2,j3; //k=i*Nv*Nv*Nv + (j1*Nv*Nv + j2*Nv + j3);
	double dx_val;
	double tmp1=0., tmp2=0., tmp3=0.;
	a[0]=0.; a[1]=0.; a[2]=0.; // the three momentum
	#pragma omp parallel for private(i,k,k0,k_v,j1,j2,j3,dx_val) shared(U) reduction(+:tmp1, tmp2, tmp3)  //reduction directive may change the result a little bit
	for(i=0; i<Nx; i++)
	{
		dx_val = dx_value(i);
		k0 = i*size_v;

		for(k_v=0; k_v<size_v; k_v++)
		{
			k = k0 + k_v;
			j3=k_v%Nv; j2=((k_v-j3)%(Nv*Nv))/Nv; j1=(k_v-j3-j2*Nv)/(Nv*Nv);
			tmp1 += Gridv((double)j1)*dx_val*U[k*6+0] + U[k*6+2]*dx_val*dv/12. + U[k*6+5]*Gridv((double)j1)*dx_val/4.;
			tmp2 += Gridv((double)j2)*dx_val*U[k*6+0] + U[k*6+3]*dx_val*dv/12. + U[k*6+5]*Gridv((double)j2)*dx_val/4.;
			tmp3 += Gridv((double)j3)*dx_val*U[k*6+0] + U[k*6+4]*dx_val*dv/12. + U[k*6+5]*Gridv((double)j3)*dx_val/4.;
		}
	}
	a[0]=tmp1*scalev; a[1]=tmp2*scalev; a[2]=tmp3*scalev;
}

void computeMomentum_Inhomo_Uniform(double *U, double *a)
{
  int k, i, j,j1,j2,j3; //k=i*Nv*Nv*Nv + (j1*Nv*Nv + j2*Nv + j3);
  double tmp1=0., tmp2=0., tmp3=0.;
  a[0]=0.; a[1]=0.; a[2]=0.; // the three momentum
  #pragma omp parallel for private(k,i,j,j1,j2,j3) shared(U) reduction(+:tmp1, tmp2, tmp3)  //reduction directive may change the result a little bit
  for(k=0;k<Nx*size_v;k++){
    j=k%size_v; i=(k-j)/size_v;
    j3=j%Nv; j2=((j-j3)%(Nv*Nv))/Nv; j1=(j-j3-j2*Nv)/(Nv*Nv);
    tmp1 += Gridv((double)j1)*dv*U[k*6+0] + U[k*6+2]*dv*dv/12. + U[k*6+5]*Gridv((double)j1)*dv/4.;
    tmp2 += Gridv((double)j2)*dv*U[k*6+0] + U[k*6+3]*dv*dv/12. + U[k*6+5]*Gridv((double)j2)*dv/4.;
    tmp3 += Gridv((double)j3)*dv*U[k*6+0] + U[k*6+4]*dv*dv/12. + U[k*6+5]*Gridv((double)j3)*dv/4.;
  }
  a[0]=tmp1*dx*dv*dv; a[1]=tmp2*dx*dv*dv; a[2]=tmp3*dx*dv*dv;
}

void computeMomentum_Homo(double *U, double *a)
{
  int k, j,j1,j2,j3; //k=i*Nv*Nv*Nv + (j1*Nv*Nv + j2*Nv + j3);
  double tmp1=0., tmp2=0., tmp3=0.;
  a[0]=0.; a[1]=0.; a[2]=0.; // the three momentum
  #pragma omp parallel for private(k,j,j1,j2,j3) shared(U) reduction(+:tmp1, tmp2, tmp3)  //reduction directive may change the result a little bit
  for(k=0;k<size_v;k++){
    j=k%size_v;
    j3=j%Nv; j2=((j-j3)%(Nv*Nv))/Nv; j1=(j-j3-j2*Nv)/(Nv*Nv);
    tmp1 += Gridv((double)j1)*dv*U[k*6+0] + U[k*6+2]*dv*dv/12. + U[k*6+5]*Gridv((double)j1)*dv/4.;
    tmp2 += Gridv((double)j2)*dv*U[k*6+0] + U[k*6+3]*dv*dv/12. + U[k*6+5]*Gridv((double)j2)*dv/4.;
    tmp3 += Gridv((double)j3)*dv*U[k*6+0] + U[k*6+4]*dv*dv/12. + U[k*6+5]*Gridv((double)j3)*dv/4.;
  }
  a[0]=tmp1*dv*dv; a[1]=tmp2*dv*dv; a[2]=tmp3*dv*dv;
}

double computeBulkVelocity_v1_in_x(double *U, int i, double x)
{
	int k, k0, k_v, j1, j2, j3;
	double dx_val, phi1_val, rho_val;
	double tmp=0.;

	dx_val = dx_value(i);
	phi1_val = (x -  Gridx((double)i))/dx_val;

	k0 = i*size_v;

	#pragma omp parallel for private(k_v,k,j1,j2,j3) shared(U) reduction(+:tmp)  //reduction directive may change the result a little bit
	for(k_v=0;k_v<size_v;k_v++)
	{
		k = k0 + k_v;
		j3=k_v%Nv; j2=((k_v-j3)%(Nv*Nv))/Nv; j1=(k_v-j3-j2*Nv)/(Nv*Nv);
		tmp += Gridv((double)j1)*(U[k*6+0] + U[k*6+1]*phi1_val) + U[k*6+2]*dv/12. + U[k*6+5]*Gridv((double)j1)/4.;
	}

	rho_val = computeMass_in_x(U, i, x);

	return tmp*scalev/(3.*rho_val);
}

double computeBulkVelocity_v2_in_x(double *U, int i, double x)
{
	int k, k0, k_v, j1, j2, j3;
	double dx_val, phi1_val, rho_val;
	double tmp=0.;

	dx_val = dx_value(i);
	phi1_val = (x -  Gridx((double)i))/dx_val;

	k0 = i*size_v;

	#pragma omp parallel for private(k_v,k,j1,j2,j3) shared(U) reduction(+:tmp)  //reduction directive may change the result a little bit
	for(k_v=0;k_v<size_v;k_v++)
	{
		k = k0 + k_v;
		j3=k_v%Nv; j2=((k_v-j3)%(Nv*Nv))/Nv; j1=(k_v-j3-j2*Nv)/(Nv*Nv);
		tmp += Gridv((double)j2)*(U[k*6+0] + U[k*6+1]*phi1_val) + U[k*6+3]*dv/12. + U[k*6+5]*Gridv((double)j2)/4.;
	}

	rho_val = computeMass_in_x(U, i, x);

	return tmp*scalev/(3.*rho_val);
}

double computeBulkVelocity_v3_in_x(double *U, int i, double x)
{
	int k, k0, k_v, j1, j2, j3;
	double dx_val, phi1_val, rho_val;
	double tmp=0.;

	dx_val = dx_value(i);
	phi1_val = (x -  Gridx((double)i))/dx_val;

	k0 = i*size_v;

	#pragma omp parallel for private(k_v,k,j1,j2,j3) shared(U) reduction(+:tmp)  //reduction directive may change the result a little bit
	for(k_v=0;k_v<size_v;k_v++)
	{
		k = k0 + k_v;
		j3=k_v%Nv; j2=((k_v-j3)%(Nv*Nv))/Nv; j1=(k_v-j3-j2*Nv)/(Nv*Nv);
		tmp += Gridv((double)j3)*(U[k*6+0] + U[k*6+1]*phi1_val) + U[k*6+4]*dv/12. + U[k*6+5]*Gridv((double)j3)/4.;
	}

	rho_val = computeMass_in_x(U, i, x);

	return tmp*scalev/(3.*rho_val);
}

double computeBulkMomentum_nu1_in_x(double *U, int i, double x)
{
	int k, k0, k_v, j1, j2, j3;
	double dx_val, phi1_val, rho_val;
	double tmp=0.;

	dx_val = dx_value(i);
	phi1_val = (x -  Gridx((double)i))/dx_val;

	k0 = i*size_v;

	#pragma omp parallel for private(k_v,k,j1,j2,j3) shared(U) reduction(+:tmp)  //reduction directive may change the result a little bit
	for(k_v=0;k_v<size_v;k_v++)
	{
		k = k0 + k_v;
		j3=k_v%Nv; j2=((k_v-j3)%(Nv*Nv))/Nv; j1=(k_v-j3-j2*Nv)/(Nv*Nv);
		tmp += Gridv((double)j1)*(U[k*6+0] + U[k*6+1]*phi1_val) + U[k*6+2]*dv/12. + U[k*6+5]*Gridv((double)j1)/4.;
	}

//	rho_val = computeMass_in_x(U, i, x);

	return tmp*scalev;	// NOTE: Was /3. (check calculation)
}

double computeBulkMomentum_nu2_in_x(double *U, int i, double x)
{
	int k, k0, k_v, j1, j2, j3;
	double dx_val, phi1_val, rho_val;
	double tmp=0.;

	dx_val = dx_value(i);
	phi1_val = (x -  Gridx((double)i))/dx_val;

	k0 = i*size_v;

	#pragma omp parallel for private(k_v,k,j1,j2,j3) shared(U) reduction(+:tmp)  //reduction directive may change the result a little bit
	for(k_v=0;k_v<size_v;k_v++)
	{
		k = k0 + k_v;
		j3=k_v%Nv; j2=((k_v-j3)%(Nv*Nv))/Nv; j1=(k_v-j3-j2*Nv)/(Nv*Nv);
		tmp += Gridv((double)j2)*(U[k*6+0] + U[k*6+1]*phi1_val) + U[k*6+3]*dv/12. + U[k*6+5]*Gridv((double)j2)/4.;
	}

//	rho_val = computeMass_in_x(U, i, x);

	return tmp*scalev;	// NOTE: Was /3. (check calculation)
}

double computeBulkMomentum_nu3_in_x(double *U, int i, double x)
{
	int k, k0, k_v, j1, j2, j3;
	double dx_val, phi1_val, rho_val;
	double tmp=0.;

	dx_val = dx_value(i);
	phi1_val = (x -  Gridx((double)i))/dx_val;

	k0 = i*size_v;

	#pragma omp parallel for private(k_v,k,j1,j2,j3) shared(U) reduction(+:tmp)  //reduction directive may change the result a little bit
	for(k_v=0;k_v<size_v;k_v++)
	{
		k = k0 + k_v;
		j3=k_v%Nv; j2=((k_v-j3)%(Nv*Nv))/Nv; j1=(k_v-j3-j2*Nv)/(Nv*Nv);
		tmp += Gridv((double)j3)*(U[k*6+0] + U[k*6+1]*phi1_val) + U[k*6+4]*dv/12. + U[k*6+5]*Gridv((double)j3)/4.;
	}

//	rho_val = computeMass_in_x(U, i, x);

	return tmp*scalev;	// NOTE: Was /3. (check calculation)
}

double computeKiE(double *U)
{
	if(Homogeneous)
	{
		return computeKiE_Homo(U);
	}
	else
	{
		return computeKiE_Inhomo(U);
	}
}

double computeKiE_Inhomo(double *U)
{
	if(MeshRefinement)
	{
		return computeKiE_Inhomo_Refined(U);
	}
	else
	{
		return computeKiE_Inhomo_Uniform(U);
	}
}

double computeKiE_Inhomo_Uniform(double *U)
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

double computeKiE_Inhomo_Refined(double *U)
{
	int k, k0, k_v, i, j,j1,j2,j3;
	double dx_val;
	double tmp=0., tp=0., tp1=0.;
	//#pragma omp parallel for private(k,i,j,j1,j2,j3,tp, tp1) shared(U) reduction(+:tmp)
	for(i=0; i<Nx; i++)
	{
		dx_val = dx_value(i);
		k0 = i*size_v;

		for(k_v=0; k_v<size_v; k_v++)
		{
			k = k0 + k_v;

			j3=k_v%Nv; j2=((k_v-j3)%(Nv*Nv))/Nv; j1=(k_v-j3-j2*Nv)/(Nv*Nv);
			//tp = ( pow(Gridv(j1+0.5), 3)- pow(Gridv(j1-0.5), 3) + pow(Gridv(j2+0.5), 3)- pow(Gridv(j2-0.5), 3) + pow(Gridv(j3+0.5), 3)- pow(Gridv(j3-0.5), 3) )/3.;
			tp1 = Gridv((double)j1)*Gridv((double)j1) + Gridv((double)j2)*Gridv((double)j2) + Gridv((double)j3)*Gridv((double)j3);
			//tmp += U[k][0]*tp + (Gridv(j1)*U[k][2]+Gridv(j2)*U[k][3]+Gridv(j3)*U[k][4])*dv*dv/6. + U[k][5]*( dv*(dv*dv*3./80. + tp1/12.) + tp/6. );
			tmp += U[k*6+0]*(tp1 + dv*dv/4.)*dx_val + (Gridv(j1)*U[k*6+2]+Gridv(j2)*U[k*6+3]+Gridv(j3)*U[k*6+4])*dx_val*dv/6. + U[k*6+5]*dx_val*(dv*dv*19./240. + tp1/4.);
		}
	}
	tmp *= scalev;
	return 0.5*tmp;
}

double computeKiE_Homo(double *U)			// int f |v|^2 dv
{
  int k, j,j1,j2,j3;
  double tmp=0., tp=0., tp1=0.;
  //#pragma omp parallel for private(k,j,j1,j2,j3,tp, tp1) shared(U) reduction(+:tmp)
  for(k=0;k<size_v;k++){
    j=k%size_v;
    j3=j%Nv; j2=((j-j3)%(Nv*Nv))/Nv; j1=(j-j3-j2*Nv)/(Nv*Nv);
    //tp = ( pow(Gridv(j1+0.5), 3)- pow(Gridv(j1-0.5), 3) + pow(Gridv(j2+0.5), 3)- pow(Gridv(j2-0.5), 3) + pow(Gridv(j3+0.5), 3)- pow(Gridv(j3-0.5), 3) )/3.;
    tp1 = Gridv((double)j1)*Gridv((double)j1) + Gridv((double)j2)*Gridv((double)j2) + Gridv((double)j3)*Gridv((double)j3);
    //tmp += U[k][0]*tp + (Gridv(j1)*U[k][2]+Gridv(j2)*U[k][3]+Gridv(j3)*U[k][4])*dv*dv/6. + U[k][5]*( dv*(dv*dv*3./80. + tp1/12.) + tp/6. );
    tmp += U[k*6+0]*(tp1 + dv*dv/4.)*dv + (Gridv(j1)*U[k*6+2]+Gridv(j2)*U[k*6+3]+Gridv(j3)*U[k*6+4])*dv*dv/6. + U[k*6+5]*( dv*dv*dv*19./240. + tp1*dv/4.);
  }
  tmp *= dv*dv;
  return 0.5*tmp;
}

double computeT_in_x(double *U, int i, double x)
{
	int k, k0, k_v, j1, j2, j3;
	double dx_val, phi1_val, rho_val, V1_val, V2_val, V3_val, V_abs_val;
	double tmp=0., tp=0., tp1=0.;

	dx_val = dx_value(i);
	phi1_val = (x -  Gridx((double)i))/dx_val;

	k0 = i*Nv*Nv*Nv;
  //#pragma omp parallel for private(k,i,j,j1,j2,j3,tp, tp1) shared(U) reduction(+:tmp)
	for(k_v=0;k_v<size_v;k_v++)
	{
		k = k0 + k_v;
		j3=k_v%Nv; j2=((k_v-j3)%(Nv*Nv))/Nv; j1=(k_v-j3-j2*Nv)/(Nv*Nv);
		//tp = ( pow(Gridv(j1+0.5), 3)- pow(Gridv(j1-0.5), 3) + pow(Gridv(j2+0.5), 3)- pow(Gridv(j2-0.5), 3) + pow(Gridv(j3+0.5), 3)- pow(Gridv(j3-0.5), 3) )/3.;
		tp1 = Gridv((double)j1)*Gridv((double)j1) + Gridv((double)j2)*Gridv((double)j2) + Gridv((double)j3)*Gridv((double)j3);
		//tmp += U[k][0]*tp + (Gridv(j1)*U[k][2]+Gridv(j2)*U[k][3]+Gridv(j3)*U[k][4])*dv*dv/6. + U[k][5]*( dv*(dv*dv*3./80. + tp1/12.) + tp/6. );
		tmp += (U[k*6+0] + U[k*6+1]*phi1_val)*(tp1 + dv*dv/4.) + (Gridv(j1)*U[k*6+2]+Gridv(j2)*U[k*6+3]+Gridv(j3)*U[k*6+4])*dv/6. + U[k*6+5]*( dv*dv*19./240. + tp1/4.);
	}

	rho_val = computeMass_in_x(U, i, x);
	V1_val = computeBulkVelocity_v1_in_x(U, i, x);
	V2_val = computeBulkVelocity_v2_in_x(U, i, x);
	V3_val = computeBulkVelocity_v3_in_x(U, i, x);
	V_abs_val = V1_val*V1_val + V2_val*V2_val + V3_val*V3_val;

	return tmp*scalev/(3.*rho_val) - V_abs_val/3.;
}

double computeKiE_in_x(double *U, int i, double x)
{
	int k, k0, k_v, j1, j2, j3;
	double dx_val, phi1_val, rho_val, V1_val, V2_val, V3_val, V_abs_val;
	double tmp=0., tp=0., tp1=0.;

	dx_val = dx_value(i);
	phi1_val = (x -  Gridx((double)i))/dx_val;

	k0 = i*Nv*Nv*Nv;
  //#pragma omp parallel for private(k,i,j,j1,j2,j3,tp, tp1) shared(U) reduction(+:tmp)
	for(k_v=0;k_v<size_v;k_v++)
	{
		k = k0 + k_v;
		j3=k_v%Nv; j2=((k_v-j3)%(Nv*Nv))/Nv; j1=(k_v-j3-j2*Nv)/(Nv*Nv);
		//tp = ( pow(Gridv(j1+0.5), 3)- pow(Gridv(j1-0.5), 3) + pow(Gridv(j2+0.5), 3)- pow(Gridv(j2-0.5), 3) + pow(Gridv(j3+0.5), 3)- pow(Gridv(j3-0.5), 3) )/3.;
		tp1 = Gridv((double)j1)*Gridv((double)j1) + Gridv((double)j2)*Gridv((double)j2) + Gridv((double)j3)*Gridv((double)j3);
		//tmp += U[k][0]*tp + (Gridv(j1)*U[k][2]+Gridv(j2)*U[k][3]+Gridv(j3)*U[k][4])*dv*dv/6. + U[k][5]*( dv*(dv*dv*3./80. + tp1/12.) + tp/6. );
		tmp += (U[k*6+0] + U[k*6+1]*phi1_val)*(tp1 + dv*dv/4.) + (Gridv(j1)*U[k*6+2]+Gridv(j2)*U[k*6+3]+Gridv(j3)*U[k*6+4])*dv/6. + U[k*6+5]*( dv*dv*19./240. + tp1/4.);
	}

	rho_val = computeMass_in_x(U, i, x);
	V1_val = computeBulkVelocity_v1_in_x(U, i, x);
	V2_val = computeBulkVelocity_v2_in_x(U, i, x);
	V3_val = computeBulkVelocity_v3_in_x(U, i, x);
	V_abs_val = V1_val*V1_val + V2_val*V2_val + V3_val*V3_val;

	return tmp*scalev/3. - V_abs_val*rho_val/3.;
}

double computeKiEratio(double *U, int *NegVals)
{
	if(Homogeneous)
	{
		return computeKiEratio_Homo(U, NegVals);
	}
	else
	{
		return computeKiEratio_Inhomo(U, NegVals);
	}
}

double computeKiEratio_Inhomo(double *U, int *NegVals)
{
	int k, k0, k_v, i,j1,j2,j3;
	double dx_val, KiEpos, KiEneg;																	// declare KiEpos (the kinetic energy where f is positive all over the cells) & KiEneg (the kinetic energy where f goes negative in the cells)
	double tmp=0., tp=0., tp1=0.;
	KiEpos = 0;
	KiEneg = 0;
	#pragma omp parallel for private(i,k,k0,k_v,j1,j2,j3,tp,tp1,dx_val) shared(U) reduction(+:KiEpos,KiEneg)
	for(i=0; i<Nx; i++)
	{
		dx_val = dx_value(i);
		k0 = i*size_v;

		for(k_v=0; k_v<size_v; k_v++)
		{
			k = k0 + k_v;
			if(NegVals[k] == 0)
			{
				j3=k_v%Nv; j2=((k_v-j3)%(Nv*Nv))/Nv; j1=(k_v-j3-j2*Nv)/(Nv*Nv);
				tp1 = Gridv((double)j1)*Gridv((double)j1) + Gridv((double)j2)*Gridv((double)j2) + Gridv((double)j3)*Gridv((double)j3);
				KiEpos += U[k*6+0]*(tp1 + dv*dv/4.)*dx_val + (Gridv(j1)*U[k*6+2]+Gridv(j2)*U[k*6+3]+Gridv(j3)*U[k*6+4])*dv*dx_val/6. + U[k*6+5]*dx_val*(dv*dv*19./240. + tp1/4.);
			}
			else
			{
				j3=k_v%Nv; j2=((k_v-j3)%(Nv*Nv))/Nv; j1=(k_v-j3-j2*Nv)/(Nv*Nv);
				tp1 = Gridv((double)j1)*Gridv((double)j1) + Gridv((double)j2)*Gridv((double)j2) + Gridv((double)j3)*Gridv((double)j3);
				KiEneg += U[k*6+0]*(tp1 + dv*dv/4.)*dx_val + (Gridv(j1)*U[k*6+2]+Gridv(j2)*U[k*6+3]+Gridv(j3)*U[k*6+4])*dv*dx_val/6. + U[k*6+5]*dx_val*(dv*dv*19./240. + tp1/4.);
			}
		}
	}
	return KiEneg/KiEpos;
}

double computeKiEratio_Homo(double *U, int *NegVals)
{
	int k, j,j1,j2,j3;
	double KiEpos, KiEneg;																	// declare KiEpos (the kinetic energy where f is positive all over the cells) & KiEneg (the kinetic energy where f goes negative in the cells)
	double tmp=0., tp=0., tp1=0.;
	KiEpos = 0;
	KiEneg = 0;
	#pragma omp parallel for private(k,j,j1,j2,j3,tp, tp1) shared(U) reduction(+:KiEpos,KiEneg)
	for(k=0;k<size_v;k++)
	{
		if(NegVals[k] == 0)
		{
			j=k%size_v;
			j3=j%Nv; j2=((j-j3)%(Nv*Nv))/Nv; j1=(j-j3-j2*Nv)/(Nv*Nv);
			tp1 = Gridv((double)j1)*Gridv((double)j1) + Gridv((double)j2)*Gridv((double)j2) + Gridv((double)j3)*Gridv((double)j3);
			KiEpos += U[k*6+0]*(tp1 + dv*dv/4.)*dv + (Gridv(j1)*U[k*6+2]+Gridv(j2)*U[k*6+3]+Gridv(j3)*U[k*6+4])*dv*dv/6. + U[k*6+5]*( dv*dv*dv*19./240. + tp1*dv/4.);
		}
		else
		{
			j=k%size_v;
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
  double ce1, cp1, dx_val;
  ce1 = computePhi_x_0(U);

  tmp1 = ce1*ce1*Lx;
  tmp2 = Lx*Lx*Lx/3.; tmp3 = -ce1*Lx*Lx;

  //#pragma omp parallel for private(j,k, i, tp1, tp2, cp1, c) shared(U) reduction(+:tmp4, tmp5, tmp6)
  for(i=0;i<Nx;i++){
	dx_val = dx_value(i);
    c = Int_Int_rho(U,i);
    cp1 = computeC_rho(U,i);
    tmp4 += dx_val*cp1 + c;
    tmp5 += dx_val*Gridx((double)i)*cp1;
    tp1=0.; tp2=0.;
    for(j=0;j<size_v;j++){
      k = i*size_v + j;
      tp1 += (U[k*6+0] + U[k*6+5]/4.);
      tp2 += U[k*6+1];
    }
    tmp5 += scalev* (tp1*( (pow(Gridx(i+0.5), 3) - pow(Gridx(i-0.5), 3))/3. - Gridx(i-0.5)*Gridx((double)i)*dx_val ) - tp2 * dx_val*dx_val*Gridx((double)i)/12.);

    tp2 *= dx_val/2.;
    tmp6 +=  cp1*cp1*dx_val + 2*cp1*c + pow(dv, 6)* ( tp1*tp1*dx_val*dx_val*dx_val/3. + tp2*tp2*dx_val/30. - tp1*tp2*dx_val*dx_val/6.);//+ tp1*tp2*(dx_val*Gridx((double)i)/6. - dx_val*dx_val/4.) ); //Int_Cumulativerho_sqr(i);
  }
  retn = tmp1 + tmp2 + tmp3 + 2*ce1*tmp4 - 2*tmp5 + tmp6;
  return 0.5*retn;
}
