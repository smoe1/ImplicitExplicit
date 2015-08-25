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

    const int numpts=xpts.getsize();

    // Water height and velocity
    double h,u;

    // Initial conditions
    for(int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);

        if (x<0.5)
        {
            h = 3.0;
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
