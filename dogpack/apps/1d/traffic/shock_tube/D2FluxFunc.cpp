#include "tensors.h"

// This is a user-supplied routine that sets the
// Derivative of the Jacobian of the Flux Function.  
// In 1d this is a 3-tensor
//
// In 1d this is a 3-tensor of size numpts x meqn x meqn x meqn
//
//     D2flux(:, i, j, k) = \partial^2 f_i / \partial q_j \partial q_k
//
// Simple Traffic Model: f(q) = q - q^2, so f''(q) = -2.0.
//
void D2FluxFunc(const dTensor1& xpts, 
        const dTensor2& Q,
        const dTensor2& Aux, 
        dTensor4& D2flux)
{
    const int numpts=xpts.getsize();
    for(int i=1; i<=numpts; i++)
    {
        // 2nd Derivative of Flux function
        D2flux.set(i, 1, 1, 1, -2.0 );

    }    
}
