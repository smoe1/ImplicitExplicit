#include "dogdefs.h"
#include "math.h"
#include "InitialParams.h"

// This is a user-supplied routine that sets the
// auxiliary arrays at all the points "xpts"
//
void AuxFunc(const dTensor2& xpts, 
	     dTensor2& auxvals)
{
  const int    numpts = xpts.getsize(1);
  const double      u = initialParams.u;
  const double      v = initialParams.v;
  const double      w = initialParams.w;

  for (int i=1; i<=numpts; i++)
    {
      double x = xpts.get(i,1);
      double y = xpts.get(i,2);
      double z = xpts.get(i,3);

      // u:  1-component of the advection velocity
      auxvals.set(i,1, u );
      
      // v:  2-component of the advection velocity
      auxvals.set(i,2, v );
      
      // w:  3-component of the advection velocity
      auxvals.set(i,3, w );      
    }

}
