#include "tensors.h"
#include "QuadMomentParams.h"


// This is a user-supplied routine that sets the source term
void SourceTermFunc(const dTensor1& xpts, 
		    const dTensor2& qvals, 
		    const dTensor2& auxvals,
                    dTensor2& fvals)
{
  const int numpts=xpts.getsize();

  // get collision time
  const double epsilon = quadMomentParams.epsilon;

  
  for (int i=1; i<=numpts; i++)
    {
      double x = xpts.get(i);

      for (int k=1; k<=fvals.getsize(2); k++)
         {
            fvals.set(i,k, 0.0 );
         }
    }



}
