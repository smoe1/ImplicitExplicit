#include<cmath>
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

    const int numpts=xpts.getsize();
    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i);
        qvals.set(i,1, exp(-100.0e0*pow(x-0.5e0,2)) * sin(80.0e0*x) );
    }
}
