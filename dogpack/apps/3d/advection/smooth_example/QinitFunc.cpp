#include <cmath>
#include "dogdefs.h"
#include "InitialParams.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, 
	       dTensor2& qvals)
{
  const int    numpts = xpts.getsize(1);
  const double     x0 = initialParams.x0;
  const double     y0 = initialParams.y0;
  const double     z0 = initialParams.z0;
  const double  width = initialParams.width;
  const double half_pi_owidth = 0.5*pi/width;

  for (int i=1; i<=numpts; i++)
    {
      const double x = xpts.get(i,1);
      const double y = xpts.get(i,2);
      const double z = xpts.get(i,3);
      
      const double r2 = pow(x-x0,2)+pow(y-y0,2)+pow(z-z0,2);
      const double r  = sqrt(r2);

      if (r<width)
        {  qvals.set(i,1, pow( cos(half_pi_owidth*r) ,6) );  }
      else
        {  qvals.set(i,1, 0.0 ); }
    }
}
