#include "dogdefs.h"
#include <cmath>
#include "QuadMomentParams.h"

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

    // Initial conditions
    for( int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i);
        double rho,u,p,q,r;

        if(x<0.5e0)
        {
            rho =  1.5;
            u   = -0.5;
            p   =  1.5;
            q   =  1.0;
            r   =  4.5;
        }
        else 
        { 
            rho =  1.0;
            u   = -0.5;
            p   =  1.0;
            q   =  0.5;
            r   =  3.0;
        }

        const double u2 = u*u;
        const double u3 = u2*u;
        const double u4 = u3*u;

        const double M2 = p + rho*u2;
        const double M3 = q + 3.0*p*u + rho*u3;
        const double M4 = r + 4.0*q*u + 6.0*p*u2 + rho*u4;

        qvals.set(i,1, rho   );  // density
        qvals.set(i,2, rho*u );  // momentum
        qvals.set(i,3, M2    );  // energy
        qvals.set(i,4, M3    );  // 3rd-order moment
        qvals.set(i,5, M4    );  // 4th-order moment
    }

}
