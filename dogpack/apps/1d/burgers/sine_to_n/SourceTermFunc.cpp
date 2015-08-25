#include "tensors.h"

// User-supplied routine that sets the source term.
//
// This is the default do-nothing source term any 1D problem.
//
void SourceTermFunc(const dTensor1& xpts, 
		    const dTensor2& qvals, 
		    const dTensor2& auxvals,
                    dTensor2& fvals)
{

    const int numpts=xpts.getsize();
    const int meqn=qvals.getsize(2);
    for(int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);
        for(int me=1; me<=meqn; me++)
        {
          fvals.set(i,me, 0.0 );
        }
    }           
        
}
