#include "dogdefs.h"

// This is a user-supplied routine that sets the
// auxiliary arrays at all the points "xpts"
//
void AuxFunc(const dTensor1& xpts, 
	     dTensor2& auxvals)
{
  const int numpts = xpts.getsize();
  const int maux   = auxvals.getsize(2);
}
