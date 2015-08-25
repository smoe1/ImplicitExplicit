#include "dogdefs.h"

// This is a user-supplied routine that sets the
// Derivative of the Jacobian of the Flux Function.  
//
// In 1d this is a 3-tensor of size numpts x meqn x meqn x meqn
//
//     D2flux(:, i, j, k) = \partial^2 f_i / \partial q_j \partial q_k
//
//     Acoustics equation
//
//           q1_t + q1_x + 3*q2_x = 0
//
//           q2_t + 3*q1_x + q2_x = 0
//
void D2FluxFunc(const dTensor1& xpts, 
		const dTensor2& Q,
		const dTensor2& Aux, 
		dTensor4& D2flux)
{
  int i,m1,m2,m3;
  int numpts=xpts.getsize();
  int meqn=Q.getsize(2);
  
  double x;
  
  for (i=1; i<=numpts; i++)
    for (m1=1; m1<=meqn; m1++)
      for (m2=1; m2<=meqn; m2++)
	for (m3=1; m3<=meqn; m3++)
	  {
	    x = xpts.get(i);
	    D2flux.set(i, m1, m2, m3, 0.0);
	  }
  
}
