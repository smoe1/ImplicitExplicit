#include "tensors.h"

// This is a user-supplied routine that sets the
// auxiliary arrays at all the points "xpts"
//
void AuxFunc(const dTensor2& xpts, 
	     dTensor2& auxvals)
{
  int numpts = xpts.getsize(1);
  int   maux = auxvals.getsize(2);
  double x,y;

   for (int i=1; i<=numpts; i++)
     for (int m=1; m<=maux; m++)
       {
	 auxvals.set(i,m, 0.0 );
       }
}
