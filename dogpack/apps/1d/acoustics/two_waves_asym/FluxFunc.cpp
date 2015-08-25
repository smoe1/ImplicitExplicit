#include "dogdefs.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     Acoustics equation
//
//           q1_t + q1_x + 3*q2_x = 0
//
//           q2_t + 3*q1_x + q2_x = 0
//
void FluxFunc(const dTensor1& xpts, 
	      const dTensor2& Q,
	      const dTensor2& Aux, 
	      dTensor2& flux)
{
  int i;
  int numpts=xpts.getsize();   
  double x;
  double h,u;

  for (i=1; i<=numpts; i++)
    {
      x = xpts.get(i);
      
      // Flux function
      flux.set(i,1, Q.get(i,1) + 3.0*Q.get(i,2) );
      flux.set(i,2, 3.0*Q.get(i,1)+Q.get(i,2) );      
    }    
}
