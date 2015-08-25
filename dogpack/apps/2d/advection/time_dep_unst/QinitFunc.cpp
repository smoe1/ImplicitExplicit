#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{

    const int numpts=xpts.getsize(1);
    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);

        const double x0 = 0.0;
        const double y0 = 0.0;

        const double r2 = pow(x-x0,2)+pow(y-y0,2);
        const double r  = sqrt(r2);

        double q = ( 0.5 + 0.5*cos( pi*r ) );
        qvals.set(i,1,q);

//      if (r<0.2)
//      {  qvals.set(i,1, pow(cos(5.0/2.0*pi*r),6) );  }
//      else
//      {  qvals.set(i,1, 0.0 ); }

        // Time starts at zero
        qvals.set(i,2,0.);
    }
}
