#include "tensors.h"
#include "TwoFluidParams.h"
#include <cmath>

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     Two-fluid plasma equations
//
void FluxFunc(const dTensor1& xpts, 
	      const dTensor2& Q,
	      const dTensor2& Aux, 
	      dTensor2& flux)
{
  const int numpts=xpts.getsize();
  const double gamma = twofluidParams.gamma;
  const double mass_ratio = twofluidParams.mass_ratio;
  const double light_speed = twofluidParams.light_speed;  
  const double ls2 = pow(light_speed,2);

  for (int i=1; i<=numpts; i++)
    {
      double x = xpts.get(i);
            
      double rho_i    = Q.get(i,1);
      double u1_i     = Q.get(i,2)/rho_i;
      double u2_i     = Q.get(i,3)/rho_i;
      double u3_i     = Q.get(i,4)/rho_i;
      double energy_i = Q.get(i,5);
      
      double rho_e    = Q.get(i,6);
      double u1_e     = Q.get(i,7)/rho_e;
      double u2_e     = Q.get(i,8)/rho_e;
      double u3_e     = Q.get(i,9)/rho_e;
      double energy_e = Q.get(i,10);
      
      double B1       = Q.get(i,11);
      double B2       = Q.get(i,12);
      double B3       = Q.get(i,13);
      double E1       = Q.get(i,14);
      double E2       = Q.get(i,15);
      double E3       = Q.get(i,16);
      
      double press_i  = (gamma-1.0e0)*(energy_i - 0.5e0*rho_i*(u1_i*u1_i 
							+ u2_i*u2_i + u3_i*u3_i));
      
      double press_e  = (gamma-1.0e0)*(energy_e - 0.5e0*rho_e*(u1_e*u1_e
							+ u2_e*u2_e + u3_e*u3_e));
      
      // Flux function
      flux.set(i,1,   rho_i*u1_i );
      flux.set(i,2,   rho_i*u1_i*u1_i + press_i );
      flux.set(i,3,   rho_i*u1_i*u2_i );
      flux.set(i,4,   rho_i*u1_i*u3_i );
      flux.set(i,5,   u1_i*(energy_i+press_i) ); 
      
      flux.set(i,6,   rho_e*u1_e );
      flux.set(i,7,   rho_e*u1_e*u1_e + press_e );
      flux.set(i,8,   rho_e*u1_e*u2_e );
      flux.set(i,9,   rho_e*u1_e*u3_e );
      flux.set(i,10,  u1_e*(energy_e+press_e) );
      
      flux.set(i,11,  0.0 );
      flux.set(i,12, -E3 );
      flux.set(i,13,  E2 );
      flux.set(i,14,  0.0 );
      flux.set(i,15,  ls2*B3 );
      flux.set(i,16, -ls2*B2 );
    }
  
}
