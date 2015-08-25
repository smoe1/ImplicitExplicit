#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{
  const int numpts=xpts.getsize(1);
  const double zfac = 1.0/sqrt(2.0*pi);
  
  for (int i=1; i<=numpts; i++)
    {
      double x = xpts.get(i,1);
      double v = xpts.get(i,2);
      
      double tmp = 1.0 + 0.04 * cos(0.3 * x);
      tmp *= zfac * ( 0.9 * exp( -0.5*pow(v,2) ) 
		      + 0.2 * exp( -4.0 * pow(v - 4.5,2) ) );
      qvals.set(i,1,  tmp);       
    }

}
