#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
// Note, for this problem: 
//
//      qex = ( 2 - cos( 2x - 2*pi*t ) ) * exp( -0.25*(4*v-1)^2 )
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{


    int i;
    int numpts=xpts.getsize(1);
    double x, v;
    double t = 0.0;
    double tmp = 0.0;

    for (i=1; i<=numpts; i++)
    {
        x = xpts.get(i,1);
        v = xpts.get(i,2);

        tmp = (2.0-cos(2.0*x)) * exp(-pow(2.0*(v-0.25),2));
        qvals.set(i, 1, tmp );
    }

    }
