#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     Simple advection equation
//
void FluxFunc(const dTensor2& xpts, const dTensor2& Q, 
        const dTensor2& Aux, dTensor3& flux)
{
    const int numpts=xpts.getsize(1);
    for(int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);

        // Variables
        const double qc = Q.get(i,1);
        const double u  = Aux.get(i,1);
        const double v  = Aux.get(i,2);

        // Flux function
        flux.set(i,1,1, u*qc );
        flux.set(i,1,2, v*qc );

//      flux.set(i,1,1, 0. );
//      flux.set(i,1,2, 0. );

        // equation #2 is used to store time only:
        flux.set(i,2,1, 0.);
        flux.set(i,2,2, 0.);

    }

}
