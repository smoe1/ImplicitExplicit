#include <cmath>
#include "dogdefs.h"
#include "InitialParams.h"
#include "EulerParams.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, 
	       dTensor2& qvals)
{
  const int    numpts = xpts.getsize(1);
  const double gamma  = eulerParams.gamma;

  for (int i=1; i<=numpts; i++)
    {
      const double x = xpts.get(i,1);
      const double y = xpts.get(i,2);
      const double z = xpts.get(i,3);
      
      double rad = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
      double rho,press,u1,u2,u3,energy;

      if(rad<0.3e0)
        {
          rho   =  initialParams.rhol;
          u1    =  initialParams.u1l;
	  u2    =  initialParams.u2l;
	  u3    =  initialParams.u3l;
	  press =  initialParams.pl;
        }
      else
        {
          rho   =  initialParams.rhor;
          u1    =  initialParams.u1r;
	  u2    =  initialParams.u2r;
	  u3    =  initialParams.u3r;
	  press =  initialParams.pr;
        }

      energy = press/(gamma-1.0e0)
        + 0.5e0*rho*(u1*u1+u2*u2+u3*u3);

      qvals.set(i,1, rho );       // density
      qvals.set(i,2, rho*u1 );    // 1-momentum
      qvals.set(i,3, rho*u2 );    // 2-momentum
      qvals.set(i,4, rho*u3 );    // 3-momentum
      qvals.set(i,5, energy );    // energy
    }
}
