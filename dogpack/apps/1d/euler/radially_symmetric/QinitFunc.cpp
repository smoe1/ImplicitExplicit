#include "tensors.h"
#include "EulerParams.h"
#include <cmath>

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

    // Gas constant
    const double gamma = eulerParams.gamma;

    // Initial conditions
    for (int i=1; i<=numpts; i++)
    {
        double rad = xpts.get(i);

        double rho,press,u1,u2,u3,energy;

        if(rad>0.2e0)
        {
            rho   =  4.0;
            u1    =  0.0;
            u2    =  0.0;
            u3    =  0.0;
            press =  1.0;
        }
        else
        {
            rho   =  4.0;
            u1    =  0.0;
            u2    =  0.0;
            u3    =  0.0;
            press =  4.0;
        }

        energy = press/(gamma-1.0e0) 
            + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);

        qvals.set(i,1, rho );       // density
        qvals.set(i,2, rho*u1 );    // 1-momentum
        qvals.set(i,3, rho*u2 );    // 2-momentum
        qvals.set(i,4, rho*u3 );    // 3-momentum
        qvals.set(i,5, energy );    // energy

    }
}
