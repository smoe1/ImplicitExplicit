#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// auxiliary arrays at all the points "xpts"
//
void AuxFunc(const dTensor2& xpts, 
	     dTensor2& auxvals)
{
  int i;
  int numpts=xpts.getsize(1);
  double x,y;
  
  for (i=1; i<=numpts; i++)
    {
      x = xpts.get(i,1);
      y = xpts.get(i,2);
	
      // u:  1-component of the advection velocity
      auxvals.set(i,1, cos(30.0*y)*exp(5*y) );
      auxvals.set(i,1, exp(5*y) );
//      auxvals.set(i,1, 100.0*y );
//      auxvals.set(i,1, 10.0*y );
//     auxvals.set(i,1, 1.0 );
     
      // v:  2-component of the advection velocity
      auxvals.set(i,2, 0.0 ); 

      // divu:  u_x + v_y
      auxvals.set(i,3,   0.0 );
    }
}
