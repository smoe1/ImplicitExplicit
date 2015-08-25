#include "tensors.h"
#include <iostream>

// This is a user-supplied routine that sets the source term
void SourceTermFunc(const dTensor1& xpts, 
		    const dTensor2& qvals, 
		    const dTensor2& auxvals,
                    dTensor2& source)
{
    int i;
    int numpts=xpts.getsize();
    double x;
    double gamma,mass_ratio,debye,cs_light,larmor_radius;
    double rho_i,u1_i,u2_i,u3_i;
    double rho_e,u1_e,u2_e,u3_e;
    double B1,B2,B3,E1,E2,E3;
    
    // Parameters
    gamma         = auxvals.get(1,1);
    mass_ratio    = auxvals.get(1,2);
    debye         = auxvals.get(1,3);
    cs_light      = auxvals.get(1,4);
    larmor_radius = auxvals.get(1,5);

    // Loop over each point
    for (i=1; i<=numpts; i++)
    {
	x = xpts.get(i);

	// Variables
	rho_i    = qvals.get(i,1);
	u1_i     = qvals.get(i,2)/rho_i;
	u2_i     = qvals.get(i,3)/rho_i;
	u3_i     = qvals.get(i,4)/rho_i;
	
	rho_e    = qvals.get(i,6);
	u1_e     = qvals.get(i,7)/rho_e;
	u2_e     = qvals.get(i,8)/rho_e;
	u3_e     = qvals.get(i,9)/rho_e;
	
	B1       = qvals.get(i,11);
	B2       = qvals.get(i,12);
	B3       = qvals.get(i,13);
	E1       = qvals.get(i,14);
	E2       = qvals.get(i,15);
	E3       = qvals.get(i,16);
	
	// Source term
	source.set(i,1,   0.0 );
	source.set(i,2,   rho_i*(E1 + u2_i*B3 - u3_i*B2)/larmor_radius );
	source.set(i,3,   rho_i*(E2 + u3_i*B1 - u1_i*B3)/larmor_radius );
	source.set(i,4,   rho_i*(E3 + u1_i*B2 - u2_i*B1)/larmor_radius );
	source.set(i,5,   rho_i*(u1_i*E1 + u2_i*E2 + u3_i*E3)/larmor_radius );
	source.set(i,6,   0.0 );
	source.set(i,7,  -rho_e*(E1 + u2_e*B3 - u3_e*B2)*mass_ratio/larmor_radius );
	source.set(i,8,  -rho_e*(E2 + u3_e*B1 - u1_e*B3)*mass_ratio/larmor_radius );
	source.set(i,9,  -rho_e*(E3 + u1_e*B2 - u2_e*B1)*mass_ratio/larmor_radius );
	source.set(i,10, -rho_e*(u1_e*E1 + u2_e*E2 + u3_e*E3)*mass_ratio/larmor_radius );
	source.set(i,11,  0.0 );
	source.set(i,12,  0.0 );
	source.set(i,13,  0.0 );
	source.set(i,14,  (mass_ratio*rho_e*u1_e-rho_i*u1_i)/(debye*debye*larmor_radius) );
	source.set(i,15,  (mass_ratio*rho_e*u2_e-rho_i*u2_i)/(debye*debye*larmor_radius) );
	source.set(i,16,  (mass_ratio*rho_e*u3_e-rho_i*u3_i)/(debye*debye*larmor_radius) );
    }
}
