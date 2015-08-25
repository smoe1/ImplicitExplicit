#include "tensors.h"

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
    int i;
    int numpts=xpts.getsize();
    double x;
    double gamma,mass_ratio,debye,cs_light,larmor_radius;
    double rho_i,u1_i,u2_i,u3_i,press_i,energy_i;
    double rho_e,u1_e,u2_e,u3_e,press_e,energy_e;
    double B1,B2,B3,E1,E2,E3;

    for (i=1; i<=numpts; i++)
    {
        x = xpts.get(i);
    
	// Parameters
	gamma         = Aux.get(i,1);
	mass_ratio    = Aux.get(i,2);
	debye         = Aux.get(i,3);
	cs_light      = Aux.get(i,4);
	larmor_radius = Aux.get(i,5);
    
	// Variables
	rho_i    = Q.get(i,1);
	u1_i     = Q.get(i,2)/rho_i;
	u2_i     = Q.get(i,3)/rho_i;
	u3_i     = Q.get(i,4)/rho_i;
	energy_i = Q.get(i,5);
	
	rho_e    = Q.get(i,6);
	u1_e     = Q.get(i,7)/rho_e;
	u2_e     = Q.get(i,8)/rho_e;
	u3_e     = Q.get(i,9)/rho_e;
	energy_e = Q.get(i,10);
	
	B1       = Q.get(i,11);
	B2       = Q.get(i,12);
	B3       = Q.get(i,13);
	E1       = Q.get(i,14);
	E2       = Q.get(i,15);
	E3       = Q.get(i,16);
    
	press_i  = (gamma-1.0e0)*(energy_i - 0.5e0*rho_i*(u1_i*u1_i 
			+ u2_i*u2_i + u3_i*u3_i));

	press_e  = (gamma-1.0e0)*(energy_e - 0.5e0*rho_e*(u1_e*u1_e
  	                + u2_e*u2_e + u3_e*u3_e));
    
	// Flux function
	flux.set(i,1,  rho_i*u1_i );
	flux.set(i,2,  rho_i*u1_i*u1_i + press_i );
	flux.set(i,3,  rho_i*u1_i*u2_i );
	flux.set(i,4,  rho_i*u1_i*u3_i );
	flux.set(i,5,  u1_i*(energy_i+press_i) ); 
    
	flux.set(i,6,  rho_e*u1_e );
	flux.set(i,7,  rho_e*u1_e*u1_e + press_e );
	flux.set(i,8,  rho_e*u1_e*u2_e );
	flux.set(i,9,  rho_e*u1_e*u3_e );
	flux.set(i,10, u1_e*(energy_e+press_e) );
	
	flux.set(i,11, 0.0 );
	flux.set(i,12, -E3 );
	flux.set(i,13,  E2 );
	flux.set(i,14, 0.0 );
	flux.set(i,15,  cs_light*cs_light*B3 );
	flux.set(i,16, -cs_light*cs_light*B2 );
    }

}
