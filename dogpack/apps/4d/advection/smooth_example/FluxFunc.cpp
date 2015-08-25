#include "dogdefs.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     Simple advection equation
//
void FluxFunc(const dTensor2& xpts,
	      const dTensor2& Q, 
              const dTensor2& Aux, 
	      dTensor3& flux)
{
  const int numpts=xpts.getsize(1);
  
  const int m = 1;
  for (int i=1; i<=numpts; i++)
    {
      const double x = xpts.get(i,1);
      const double y = xpts.get(i,2);
      const double z = xpts.get(i,3);
      
      // Variables
      const double qc = Q.get(i,m);
      const double u  = Aux.get(i,1);
      const double v  = Aux.get(i,2);
      const double w  = Aux.get(i,3);
      
      // Flux function
      flux.set(i,m,1, u*qc );
      flux.set(i,m,2, v*qc );
      flux.set(i,m,3, w*qc );
    }

}
