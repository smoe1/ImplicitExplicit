#include "dogdefs.h"
#include "dog_math.h"

double GetAlpha(const double rho,
                const double p,
                const double q,
                const double r)
{
  // check inputs
  const double p2 = p*p;
  const double p3 = p*p2;
  const double q2 = q*q;
  const double rmin = (p2/rho + q2/p) - 1.0e-8;
  const double rmax = (13.0*q2)/(3.0*p) + (33.0*p2)/(13.0*rho) + 1.0e-8;
  
  if (r>=((13.0*q2)/(3.0*p) + (33.0*p2)/(13.0*rho)))
    {
      return 3.0/13.0;
    }
  
  /*
  if (rho<=0.0 || p<=0.0 || r<rmin || r>rmax )
    {
      printf("  rho = %22.16e \n",rho);
      printf("    p = %22.16e \n",p);
      printf("    q = %22.16e \n",q);
      printf("    r = %22.16e \n",r);
      printf("\n");
      printf(" rmin = %22.16e \n",rmin);
      printf(" rmax = %22.16e \n",rmax);
      printf("\n");
      exit(1);
    }
  */
  
  double  Poly(double a, const double rho, const double p, const double q, const double r);
  double DPoly(double a, const double rho, const double p, const double q, const double r);
   
  double al = (p2/rho + q2/p);
  double ar = (13.0*q2)/(3.0*p) + (33.0*p2)/(13.0*rho);
  double fal = Poly(al,rho,p,q,r);
  double far = Poly(ar,rho,p,q,r);

  double a = (far*al - fal*ar)/(far-fal);

  int mfound = 0;

  double Pa = Poly(a,rho,p,q,r);
  if (fabs(Pa)<=1.0e-12)
    { mfound = 1; }

  int NumIters = 0;
  while(mfound==0)
    {
      double corr = Pa/DPoly(a,rho,p,q,r);
      a = a - corr;
      Pa = Poly(a,rho,p,q,r);

      if (fabs(Pa)<=1.0e-12 && fabs(corr)<=1.0e-12)
	{
	  mfound = 1;
	}
      else
	{
	  NumIters = NumIters + 1;
	}

      if (NumIters>1000)
	{
	  printf(" ERROR in GetAlpha, NumIters = %i\n",NumIters);
	  exit(1);
	}
    }

  return a;

}

double Poly(double a, const double rho, const double p, const double q, const double r)
{
  double p2 = p*p;
  double p3 = p2*p;
  double q2 = q*q;

  double a2 = a*a;
  double a3 = a2*a;

  return (13.0*p3*a3 - 6.0*p3*a2 + a*p*(-12.0*p2 + 5.0*r*rho) - 5.0*rho*q2);
}

double DPoly(double a, const double rho, const double p, const double q, const double r)
{
  double p2 = p*p;
  double p3 = p2*p;
  double q2 = q*q;

  double a2 = a*a;

  return (39.0*p3*a2 - 12.0*p3*a + p*(-12.0*p2 + 5.0*r*rho) );
}
