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
  const double     u1 = initialParams.u1;
  const double     u2 = initialParams.u2;
  const double     u3 = initialParams.u3;
  const double     u4 = initialParams.u4;

  for (int i=1; i<=numpts; i++)
    {
      double x = xpts.get(i,1);
      double y = xpts.get(i,2);
      double z = xpts.get(i,3);
      double w = xpts.get(i,4);

      // 1-component of the advection velocity
      auxvals.set(i,1, u1 );
      
      // 2-component of the advection velocity
      auxvals.set(i,2, u2 );
      
      // 3-component of the advection velocity
      auxvals.set(i,3, u3 );      

      // 3-component of the advection velocity
      auxvals.set(i,4, u4 );      

    }

}
