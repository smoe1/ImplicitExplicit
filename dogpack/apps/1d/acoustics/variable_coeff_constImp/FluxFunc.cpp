#include "dogdefs.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     Acoustics equation:  q_tt - q_xx = 0.
//
//     Substituting q1 = q_t and q2 = q2_x, in flux form,
//
//           q1_t - q2_x = 0
//
//           q2_t - q1_x = 0
//
void FluxFunc(const dTensor1& xpts, 
	      const dTensor2& Q,
	      const dTensor2& Aux, 
	      dTensor2& flux)
{

    const int numpts=xpts.getsize();   
    double rho;
    double c;
    double bk;
    for (int i=1; i<=numpts; i++)
    {
    rho=Aux.get(i,2);
    c=Aux.get(i,1);
    bk=c*c*rho;

        double x = xpts.get(i);

        // Flux function
        flux.set(i,1, bk*Q.get(i,2) );
        flux.set(i,2, 1.0/rho*Q.get(i,1) );

    }    

}
