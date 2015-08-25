#include <cmath>
#include "dogdefs.h"
#include "VlasovParams.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{

    const int numpts=xpts.getsize(1);

//  const double rho = vlasovParams.rho0;
//  const double T   = vlasovParams.temp;

//const double alpha = 0.01;
//const double k = 0.5;

    const double alpha = vlasovParams.alpha;
    const double k     = vlasovParams.k;

    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double v = xpts.get(i,2);

        double tmp = 1.0/sqrt(2.0*pi) * (1.0+alpha*cos(k*x)) * exp( -v*v / 2.0 );

        qvals.set(i,1,  tmp); 

    }

}
