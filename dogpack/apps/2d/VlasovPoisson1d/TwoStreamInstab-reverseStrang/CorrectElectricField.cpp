#include "dogdefs.h"
#include "constants.h"
#include "DogParamsCart2.h"
#include "dog_math.h"

//
// Correct electric field to guarantee energy conservation
//
void CorrectElectricField(const double dt,
			  const dTensorBC4* qold,
			  const dTensorBC4* qnew,
			  const dTensorBC3* Efield_old,
			  dTensorBC3* Efield)
{
  const int morder = Efield->getsize(3);
  const int mx = dogParamsCart2.get_mx();
  const int my = dogParamsCart2.get_my();
  const double dx = dogParamsCart2.get_dx();
  const double dy = dogParamsCart2.get_dy();
  const double dx2 = dx*dx;
  const double dy2 = dy*dy;
  const double dx3 = dx2*dx;
  const double dy3 = dy2*dy;
  const double dx4 = dx3*dx;
  const double dy4 = dy3*dy;
  const double xlow = dogParamsCart2.get_xlow();
  const double ylow = dogParamsCart2.get_ylow();

  double Gold_p[mx+2];
  double Gold_m[mx+2];
  double Gnew_p[mx+2];
  double Gnew_m[mx+2];

  for (int i=1; i<=mx; i++)
    { 
      double GMold = 0.0;
      double GPold = 0.0;
      double GMnew = 0.0;
      double GPnew = 0.0;

      double F[16];
      for (int j=1; j<=my; j++)
	{
	  double vj = ylow + (double(j)-0.5)*dy;
	  double vj2 = vj*vj;
	  double vj3 = vj2*vj;
	  double vj4 = vj3*vj;

	  for (int m=1; m<=15; m++)
	    {  F[m] = qold->get(i,j,1,m);  }

	  GMold = GMold + 
	    0.5*F[1]*vj3*dy+0.5*F[1]*vj*dy3+0.1*F[3]*dy4*sq3+0.5*F[3]*dy2*sq3*vj2+0.2*F[6]*dy3*vj*sq5+(1/35)*F[10]*dy4*sq7-0.5*sq3*F[2]*dy*vj3-0.5*sq3*F[2]*dy3*vj-0.2*F[8]*dy3*vj*sq3*sq5-(1/35)*F[12]*dy4*sq3*sq7-0.3*F[4]*dy4-1.5*F[4]*dy2*vj2+0.5*sq5*F[5]*dy*vj3+0.5*sq5*F[5]*dy3*vj+0.1*sq5*F[7]*dy4*sq3+0.5*sq5*F[7]*dy2*sq3*vj2+F[13]*dy3*vj-0.1*dy4*sq7*F[11]*sq3-0.5*dy2*sq7*F[11]*sq3*vj2-0.5*dy*sq7*F[9]*vj3-0.5*dy3*sq7*F[9]*vj+1.5*F[14]*dy*vj3+1.5*F[14]*dy3*vj;

	  GPold = GPold + 
	    0.5*F[1]*vj3*dy+0.5*F[1]*vj*dy3+0.1*F[3]*dy4*sq3+0.5*F[3]*dy2*sq3*vj2+0.2*F[6]*dy3*vj*sq5+(1/35)*F[10]*dy4*sq7+0.5*sq3*F[2]*dy*vj3+0.5*sq3*F[2]*dy3*vj+0.2*F[8]*dy3*vj*sq3*sq5+(1/35)*F[12]*dy4*sq3*sq7+0.3*F[4]*dy4+1.5*F[4]*dy2*vj2+0.5*sq5*F[5]*dy*vj3+0.5*sq5*F[5]*dy3*vj+0.1*sq5*F[7]*dy4*sq3+0.5*sq5*F[7]*dy2*sq3*vj2+F[13]*dy3*vj+0.1*dy4*sq7*F[11]*sq3+0.5*dy2*sq7*F[11]*sq3*vj2+0.5*dy*sq7*F[9]*vj3+0.5*dy3*sq7*F[9]*vj+1.5*F[14]*dy*vj3+1.5*F[14]*dy3*vj;

	  for (int m=1; m<=15; m++)
	    {  F[m] = qnew->get(i,j,1,m);  }

	  GMnew = GMnew + 
	    0.5*F[1]*vj3*dy+0.5*F[1]*vj*dy3+0.1*F[3]*dy4*sq3+0.5*F[3]*dy2*sq3*vj2+0.2*F[6]*dy3*vj*sq5+(1.0/35.0)*F[10]*dy4*sq7-0.5*sq3*F[2]*dy*vj3-0.5*sq3*F[2]*dy3*vj-0.2*F[8]*dy3*vj*sq3*sq5-(1.0/35.0)*F[12]*dy4*sq3*sq7-0.3*F[4]*dy4-1.5*F[4]*dy2*vj2+0.5*sq5*F[5]*dy*vj3+0.5*sq5*F[5]*dy3*vj+0.1*sq5*F[7]*dy4*sq3+0.5*sq5*F[7]*dy2*sq3*vj2+F[13]*dy3*vj-0.1*dy4*sq7*F[11]*sq3-0.5*dy2*sq7*F[11]*sq3*vj2-0.5*dy*sq7*F[9]*vj3-0.5*dy3*sq7*F[9]*vj+1.5*F[14]*dy*vj3+1.5*F[14]*dy3*vj;

	  GPnew = GPnew + 
	    0.5*F[1]*vj3*dy+0.5*F[1]*vj*dy3+0.1*F[3]*dy4*sq3+0.5*F[3]*dy2*sq3*vj2+0.2*F[6]*dy3*vj*sq5+(1.0/35.0)*F[10]*dy4*sq7+0.5*sq3*F[2]*dy*vj3+0.5*sq3*F[2]*dy3*vj+0.2*F[8]*dy3*vj*sq3*sq5+(1.0/35.0)*F[12]*dy4*sq3*sq7+0.3*F[4]*dy4+1.5*F[4]*dy2*vj2+0.5*sq5*F[5]*dy*vj3+0.5*sq5*F[5]*dy3*vj+0.1*sq5*F[7]*dy4*sq3+0.5*sq5*F[7]*dy2*sq3*vj2+F[13]*dy3*vj+0.1*dy4*sq7*F[11]*sq3+0.5*dy2*sq7*F[11]*sq3*vj2+0.5*dy*sq7*F[9]*vj3+0.5*dy3*sq7*F[9]*vj+1.5*F[14]*dy*vj3+1.5*F[14]*dy3*vj;
	}

      Gold_m[i] = GMold;
      Gold_p[i] = GPold;
      Gnew_m[i] = GMnew;
      Gnew_p[i] = GPnew;
    }
  Gold_p[0]    = Gold_p[mx];
  Gold_m[mx+1] = Gold_m[1];
  Gnew_p[0]    = Gnew_p[mx];
  Gnew_m[mx+1] = Gnew_m[1];

  for (int i=1; i<=mx; i++)
    {
      double Znew = 0.0;
      double Zold = 0.0;

      double Enew = 0.0;
      double Eold = 0.0;

      double Zcorr = 0.0;

      for (int j=1; j<=my; j++)
	{
	  double vj = ylow + (double(j)-0.5)*dy;
	  double vj2 = vj*vj;
	  double F1,F3,F6;

	  F1 = qnew->get(i,j,1,1);
	  F3 = qnew->get(i,j,1,3);
	  F6 = qnew->get(i,j,1,6);
	  Enew = Enew + (onethird*onehalf*dy)*(3.0*vj2+dy2)*F1
	    + (dy2*osq3)*vj*F3 + onethird*osq5*dy3*F6;

	  F1 = qold->get(i,j,1,1);
	  F3 = qold->get(i,j,1,3);
	  F6 = qold->get(i,j,1,6);
	  Eold = Eold + (onethird*onehalf*dy)*(3.0*vj2+dy2)*F1
	    + (dy2*osq3)*vj*F3 + onethird*osq5*dy3*F6;
	}
      
      for (int m=1; m<=morder; m++)
	{	  
	  Zold = Zold + 0.5*pow(Efield_old->get(i,1,m),2);
	  Znew = Znew + 0.5*pow(Efield->get(i,1,m),2);
	}

      
      Zcorr = Zold + Eold - Enew - (0.125/dx)*( (Gnew_p[i]+Gnew_m[i+1] - Gnew_p[i-1]+Gnew_m[i]) 
				             + (Gold_p[i]+Gold_m[i+1] - Gold_p[i-1]+Gold_m[i]) );

      for (int m=1; m<=morder; m++)
	{
	  Efield->set(i,1,m, sqrt(Zcorr/Znew)*Efield->get(i,1,m) );
	}
    }
}
