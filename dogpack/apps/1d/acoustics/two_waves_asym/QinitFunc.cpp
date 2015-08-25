#include <cmath>
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

    const int numpts = xpts.getsize();
    double tmp = 0.0;
    double a   = 0.2;

    // Initial conditions
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);

        if( fabs(x+0.2) < a)
        {
            tmp = pow(cos( (x+0.2) * pi/(2.0*a)),6);
            qvals.set(i, 1, pow( cos( (x+0.2) * pi/(2.0*a)), 6) );
        }
        else
        {
            qvals.set(i,1, 0.0);
        }

        qvals.set(i, 2, 0.0);

    }

}
