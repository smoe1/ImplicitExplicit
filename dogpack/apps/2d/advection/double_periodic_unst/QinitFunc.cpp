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
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        double x0 = 0.5;
        double y0 = 0.5;

        double r2 = pow(x-x0,2)+pow(y-y0,2);
        double r  = sqrt(r2);

        if (r<0.3)
        {  qvals.set(i,1, pow(cos(5.0/3.0*pi*r),6) );  }
        else
        {  qvals.set(i,1, 0.0 ); }
    }

}
