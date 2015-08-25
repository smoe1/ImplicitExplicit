#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
// Each application is REQUIRED to define one of these.
//
// Input:
//
//    xpts( 1:numpts )           - The x-coordinates for a list of points
//
// Output:
//
//    qvals( 1:numpts, 1:meqn )  - The vector of conserved variables, q at
//                                 each point.
//
// See also: AuxFunc.
void QinitFunc(const dTensor1& xpts, dTensor2& qvals)
{

  const int numpts = xpts.getsize();

  // Parameters describing initial condition.
  const double a   = 0.2;
  const double xshift = 0.0;

  // Initial conditions
  for(int i=1; i<=numpts; i++)
    {
      const double x = xpts.get(i);
      if( fabs(x+xshift) < a)
        {
	  double tmp = pow(cos( (x+xshift) * pi/(2.0*a)),6);
	  qvals.set(i, 1, pow( cos( (x+xshift) * pi/(2.0*a)), 6) );
        }
      else
        {
	  qvals.set(i,1, 0.0);
        }

      qvals.set(i, 2, 0.0);

    }

}
