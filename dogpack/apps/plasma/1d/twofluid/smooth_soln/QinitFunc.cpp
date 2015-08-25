#include "dogdefs.h"
#include "TwoFluidParams.h"
#include <cmath>

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor1& xpts, 
	       dTensor2& qvals)
{
  const int numpts = xpts.getsize();
  const double gamma = twofluidParams.gamma;
  const double mass_ratio = twofluidParams.mass_ratio;
  const double light_speed = twofluidParams.light_speed;  
    
  // Initial conditions
  for (int i=1; i<=numpts; i++)
    {
      double x = xpts.get(i);

      double nfunc = 1.0 + 0.5*sin(4.0*pi*x);
      double afunc = cos(4.0*pi*x);
      double bfunc = sin(2.0*pi*x);

      double rho_i   =  nfunc;
      double u1_i    =  1.0;
      double u2_i    =  0.0;
      double u3_i    =  0.0;
      double press_i =  1.0;
      
      double rho_e   =  nfunc/mass_ratio;
      double u1_e    =  1.0;
      double u2_e    =  0.0;
      double u3_e    =  0.0;
      double press_e =  0.2/mass_ratio;

      double B1      =  1.0;
      double B2      =  afunc;
      double B3      =  bfunc;
      double E1      =  0.0;
      double E2      =  bfunc;
      double E3      = -afunc;	    

      double energy_i = press_i/(gamma-1.0e0) 
	+ 0.5e0*rho_i*(u1_i*u1_i + u2_i*u2_i + u3_i*u3_i);
      
      double energy_e = press_e/(gamma-1.0e0) 
	+ 0.5e0*rho_e*(u1_e*u1_e + u2_e*u2_e + u3_e*u3_e);
      
      qvals.set(i,1, rho_i );       // density (ion)
      qvals.set(i,2, rho_i*u1_i );  // 1-momentum (ion)
      qvals.set(i,3, rho_i*u2_i );  // 2-momentum (ion)
      qvals.set(i,4, rho_i*u3_i );  // 3-momentum (ion)
      qvals.set(i,5, energy_i );    // energy (ion)
      
      qvals.set(i,6, rho_e );       // density (electron)
      qvals.set(i,7, rho_e*u1_e );  // 1-momentum (electron)
      qvals.set(i,8, rho_e*u2_e );  // 2-momentum (electron)
      qvals.set(i,9, rho_e*u3_e );  // 3-momentum (electron)
      qvals.set(i,10, energy_e );   // energy (electron)
      
      qvals.set(i,11, B1 );         // 1-magnetic field
      qvals.set(i,12, B2 );         // 2-magnetic field
      qvals.set(i,13, B3 );         // 3-magnetic field
      qvals.set(i,14, E1 );         // 1-electric field
      qvals.set(i,15, E2 );         // 2-electric field
      qvals.set(i,16, E3 );         // 3-electric field     
    }
}
