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
void QinitFunc(const dTensor1& xpts,
        dTensor2& qvals)
{
    int i;
    int numpts=xpts.getsize();
    double x;
    double h,u;

    double a = 0.4;

    // Initial conditions
    for (i=1; i<=numpts; i++)
    {
        x = xpts.get(i);


        if ( fabs(x-0.5) < a)
        {
            h = 1.0 + pow(cos(pi*(x-0.5) / (2.0*a) ),4);
            u = 0.0;
        }
        else
        {
            h = 1.0;
            u = 0.0;
        }

        qvals.set(i,1, h );    // depth
        qvals.set(i,2, h*u );  // momentum    
    }
}
