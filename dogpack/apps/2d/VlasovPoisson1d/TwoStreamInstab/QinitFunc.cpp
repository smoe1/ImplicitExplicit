#include <cmath>
#include "dogdefs.h"
#include "VlasovParams.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{
    const int numpts=xpts.getsize(1);

    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double v = xpts.get(i,2);

        qvals.set(i, 1, 
            pow(v,2) / sqrt(2.0*pi) * (
                1.0 - vlasovParams.alpha*cos( vlasovParams.k *x) ) * exp( -0.5 * v*v ) );
    }

}
