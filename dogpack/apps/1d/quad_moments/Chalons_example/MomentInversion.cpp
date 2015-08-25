#include "dogdefs.h"
#include "dog_math.h"

double MomentInversion(const dTensor1& Q, dTensor1& BG)
{
  const double M0 = Q.get(1);
  const double M1 = Q.get(2);
  const double M2 = Q.get(3);
  const double M3 = Q.get(4);
  const double M4 = Q.get(5);
  
  const double rho = M0;
  const double u   = M1/rho;
  const double u2  =  u*u;
  const double u3  = u2*u;
  const double u4  = u3*u;
  const double u5  = u4*u;
  const double p   = M2 - rho*u2;
  const double p2  = p*p;
  const double p3 =  p2*p;
  const double q   = M3 - 3.0*p*u - rho*u3;
  const double r   = M4 - 4.0*q*u - 6.0*p*u2 - rho*u4;
  
  double GetAlpha(const double rho,
                  const double p,
                  const double q,
                  const double r);
  const double alpha = GetAlpha(rho,p,q,r);

  double om1,om2,mu1,mu2,sig;
  if (alpha<1.0e-10)
    {
      om1 = 0.5*rho;
      om2 = 0.5*rho;
      mu1 = u;
      mu2 = u;
      sig = p/rho;
    }
  else
    {
      double qt = q/alpha;
      double qt2 = qt*qt;
      double rho2 = rho*rho;
      double tmp1 = rho2*qt/sqrt(rho2*qt2+4.0*rho*p3*alpha);

      om1 = 0.5*(rho + tmp1);
      om2 = 0.5*(rho - tmp1);

      double tmp2 = sqrt(p*alpha/rho + pow(0.5*qt/p,2));

      mu1 = u + 0.5*qt/p - tmp2;
      mu2 = u + 0.5*qt/p + tmp2;

      sig = p/rho*(1.0-alpha);
    }

  BG.set(1, om1 );
  BG.set(2, om2 );
  BG.set(3, mu1 );
  BG.set(4, mu2 );
  BG.set(5, sig );

  return alpha;
}
