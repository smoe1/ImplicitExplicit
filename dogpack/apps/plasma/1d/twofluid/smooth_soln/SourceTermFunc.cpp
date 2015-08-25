#include "tensors.h"
#include "TwoFluidParams.h"
#include <cmath>

// This is a user-supplied routine that sets the source term
void SourceTermFunc(const dTensor1& xpts, 
		    const dTensor2& qvals, 
		    const dTensor2& auxvals,
                    dTensor2& source)
{
  const int numpts=xpts.getsize();
  const double gamma = twofluidParams.gamma;
  const double mass_ratio = twofluidParams.mass_ratio;
  const double light_speed = twofluidParams.light_speed;  
  const double ls2 = pow(light_speed,2);
    
  // Loop over each point
  for (int i=1; i<=numpts; i++)
    {
      double x = xpts.get(i);
      
      // Variables
      double rho_i = qvals.get(i,1);
      double u1_i  = qvals.get(i,2)/rho_i;
      double u2_i  = qvals.get(i,3)/rho_i;
      double u3_i  = qvals.get(i,4)/rho_i;
	
      double rho_e = qvals.get(i,6);
      double u1_e  = qvals.get(i,7)/rho_e;
      double u2_e  = qvals.get(i,8)/rho_e;
      double u3_e  = qvals.get(i,9)/rho_e;
	
      double B1    = qvals.get(i,11);
      double B2    = qvals.get(i,12);
      double B3    = qvals.get(i,13);
      double E1    = qvals.get(i,14);
      double E2    = qvals.get(i,15);
      double E3    = qvals.get(i,16);
      
      // Source term
      source.set(i,1,   0.0 );
      source.set(i,2,   rho_i*(E1 + u2_i*B3 - u3_i*B2) );
      source.set(i,3,   rho_i*(E2 + u3_i*B1 - u1_i*B3) );
      source.set(i,4,   rho_i*(E3 + u1_i*B2 - u2_i*B1) );
      source.set(i,5,   rho_i*(u1_i*E1 + u2_i*E2 + u3_i*E3) );
      source.set(i,6,   0.0 );
      source.set(i,7,  -rho_e*(E1 + u2_e*B3 - u3_e*B2)*mass_ratio );
      source.set(i,8,  -rho_e*(E2 + u3_e*B1 - u1_e*B3)*mass_ratio );
      source.set(i,9,  -rho_e*(E3 + u1_e*B2 - u2_e*B1)*mass_ratio );
      source.set(i,10, -rho_e*(u1_e*E1 + u2_e*E2 + u3_e*E3)*mass_ratio );
      source.set(i,11,  0.0 );
      source.set(i,12,  0.0 );
      source.set(i,13,  0.0 );
      source.set(i,14,  ls2*(mass_ratio*rho_e*u1_e-rho_i*u1_i) );
      source.set(i,15,  ls2*(mass_ratio*rho_e*u2_e-rho_i*u2_i) );
      source.set(i,16,  ls2*(mass_ratio*rho_e*u3_e-rho_i*u3_i) );
    }
}
