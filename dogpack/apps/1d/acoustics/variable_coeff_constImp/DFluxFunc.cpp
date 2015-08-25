#include "dogdefs.h"

// This is a user-supplied routine that sets the
// Jacobian of the Flux Function 
//
//     Acoustics equation
//
void DFluxFunc(const dTensor1& xpts, 
        const dTensor2& Q,
        const dTensor2& Aux, 
        dTensor3& Dflux)
{
    int i,j;
    int numpts=xpts.getsize();
    double x;

    for (i=1; i<=numpts; i++)
    {
        x = xpts.get(i);

        Dflux.set(i,1,1, 0.0);
        Dflux.set(i,1,2, 1.0e0);
        Dflux.set(i,2,1, 1.0e0);
        Dflux.set(i,2,2, 0.0);
    }

}
