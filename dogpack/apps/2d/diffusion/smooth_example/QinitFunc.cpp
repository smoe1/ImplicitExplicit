#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{
  int i;
  int numpts=xpts.getsize(1);
  double x,y;
  
  for (i=1; i<=numpts; i++)
    {
      x = xpts.get(i,1);
      y = xpts.get(i,2);
      
      qvals.set(i,1, sin(2.0*pi*x)*sin(2.0*pi*y) );
    }
}
