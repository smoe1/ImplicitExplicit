#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(
    const dTensor2& xpts, dTensor2& qvals)
{

    const int numpts  = qvals.getsize(1);
    for (int i=1; i<=numpts; i++)
    {

        const double x  = xpts.get(i,1);
        const double y  = xpts.get(i,2);

        double r  = sqrt( pow(x-0.5,2) + pow(y-0.5,2) );

        qvals.set(i,1, 0.5 + 0.5*cos( pi*r ) );

    }
}
