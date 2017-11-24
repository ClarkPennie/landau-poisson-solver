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
