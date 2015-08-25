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
      const double r   = M4 - 4.0*q*u - 6.0*p*u2 - rho*u4;
      
      double GetAlpha(const double rho,
		      const double p,
		      const double q,
		      const double r);
      double alpha = GetAlpha(rho,p,q,r);

      double qt;
      if (alpha < 1.0e-10)
        {  qt = 0.0;  }
      else
        {  qt = q/alpha;  }
      double qt2 = qt*qt;
      double qt3 = qt2*qt;

      const double M5t = (qt3/p2) + (5.0*qt2*u/p) + 10.0*qt*u2 
        + 10.0*p*qt/rho - 2.0*p*alpha/rho*(4.0*qt+5.0*p*u);

      flux.set(i,1, M1 );
      flux.set(i,2, M2 );
      flux.set(i,3, M3 ); 
      flux.set(i,4, rho*u4 + 6.0*p*u2 + 4*q*u + alpha*qt2/p + 3.0*p2/rho - 2*p2*alpha*alpha/rho );
      flux.set(i,5, rho*u5 + 10.0*p*u3 + 15.0*p2*u/rho + alpha*M5t );
    }

}
