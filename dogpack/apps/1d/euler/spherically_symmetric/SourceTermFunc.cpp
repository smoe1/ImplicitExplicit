#include "tensors.h"
#include "EulerParams.h"

// This is a user-supplied routine that sets the source term
void SourceTermFunc(const dTensor1& xpts, 
		    const dTensor2& qvals, 
		    const dTensor2& auxvals,
                    dTensor2& fvals)
{
  const int numpts=xpts.getsize();

  // Gas constant
  const double gamma = eulerParams.gamma;

  // Geometric source term from radial symmetry
  for (int i=1; i<=numpts; i++)
    {
      double rad = xpts.get(i);
           
      // Variables
      double rho    = qvals.get(i,1);
      double u1     = qvals.get(i,2)/rho;
      double energy = qvals.get(i,3);
      double press  = (gamma-1.0e0)*(energy-0.5e0*rho*(u1*u1));
      
      // Source terms
      fvals.set(i,1,  -2.0*rho*u1/rad            );
      fvals.set(i,2,  -2.0*rho*u1*u1/rad         );
      fvals.set(i,3,  -2.0*u1*(energy+press)/rad );
    }
}
