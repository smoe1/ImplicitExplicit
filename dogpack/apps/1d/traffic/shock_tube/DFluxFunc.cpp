#include "tensors.h"

// This is a user-supplied routine that sets the
// Jacobian of the Flux Function 
// 
// The expected format is Dflux.get(:, i, j) = \partial f_i, \partial q_j.
//
//     Simple Traffic Model, f'(q) = 1 - 2q
//
void DFluxFunc(const dTensor1& xpts, 
        const dTensor2& Q,
        const dTensor2& Aux, 
        dTensor3& Dflux)
{
    const int numpts=xpts.getsize();
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);
        Dflux.set(i,1,1, 1.0 - 2.0*Q.get(i,1)  );
    }    
}
