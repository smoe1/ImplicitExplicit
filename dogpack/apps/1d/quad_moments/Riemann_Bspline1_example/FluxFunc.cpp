#include "dogdefs.h"
#include "dog_math.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
void FluxFunc(const dTensor1& xpts, 
	      const dTensor2& Qin,
	      const dTensor2& Aux, 
	      dTensor2& flux)
{
  const int numpts=xpts.getsize();
  
  for (int i=1; i<=numpts; i++)
    {
      const double M0 = Qin.get(i,1);
      const double M1 = Qin.get(i,2);
      const double M2 = Qin.get(i,3);
      const double M3 = Qin.get(i,4);
      const double M4 = Qin.get(i,5);

      const double rho = M0;
      const double u   = M1/rho;
      const double u2  =  u*u;
      const double u3  = u2*u;
      const double u4  = u3*u;
      const double u5  = u4*u;
      const double p   = M2 - rho*u2;
      const double p2  = p*p;
      const double q   = M3 - 3.0*p*u - rho*u3;      
      const double q2  = q*q;
      const double q3  = q*q2;
      const double r   = M4 - 4.0*q*u - 6.0*p*u2 - rho*u4;

      double GetAlpha(const double rho,
		      const double p,
		      const double q,
		      const double r);
            
      const double alpha = GetAlpha(rho,p,q,r);
      const double alpha2 = alpha*alpha;

      //const double M4 = rho*u4 + 6.0*p*u2 + 4.0*q*u - (13.0/5.0)*p2*alpha2/rho
      //+ (6.0/5.0)*p2*alpha/rho + (12.0/5.0)*p2/rho + q2/(p*alpha);

      const double M5bar = rho*u5 + 10.0*u2*(p*u+q) + 2.0*p*q*(5.0-4.0*alpha)/rho 
	+ p2*u*(12.0+6.0*alpha-13.0*alpha2)/rho + 5.0*q2*u/(p*alpha) + q3/(p2*alpha2);

      flux.set(i,1, M1 );
      flux.set(i,2, M2 );
      flux.set(i,3, M3 ); 
      flux.set(i,4, M4 );
      flux.set(i,5, M5bar );
    }

}
