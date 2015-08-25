#include "dogdefs.h"
// This is a user-supplied routine that sets the
// auxiliary arrays at all the points "xpts"
//
void AuxFunc(const dTensor1& xpts, 
	     dTensor2& auxvals)
{
  const int numpts=xpts.getsize();
  
  const double rhol=1.0;    
  const double cl=1.0;    
  const double rhor=2.0;    
  const double cr=0.5;  
  
  
  for (int i=1; i<=numpts; i++)
    {
      double x = xpts.get(i);
      if(x<0.0)
        {
	  auxvals.set(i,1, cl );
	  auxvals.set(i,2, rhol);
        }
      else
        {
	  auxvals.set(i,1, cr );
	  auxvals.set(i,2, rhor);
        }
    }
}
