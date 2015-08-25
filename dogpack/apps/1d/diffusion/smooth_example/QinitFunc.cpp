#include <cmath>
#include "constants.h"
#include "tensors.h"

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
    int i;
    int numpts=xpts.getsize();
    double x;

    for (i=1; i<=numpts; i++)
    {
        x = xpts.get(i);

        qvals.set(i,1, sin(4.0*pi*x) );
        //if ( fabs(x-0.25e0) > 0.125e0 )
        //  { qvals.set(i,1, 0.0e0 ); }
        //else
        //  { qvals.set(i, 1,  pow(cos(pi*(4.0e0*x-1.0e0)),6) ); }

    }

}

