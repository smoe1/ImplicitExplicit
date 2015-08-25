#include "tensors.h"
#include "Params.h"

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
    int numpts=xpts.getsize();

    // Parameters
    double cs_light_2 = appParams.cs_light*appParams.cs_light;
    
    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i);
    
	// Variables
	const double& rho_i = Q.get(i,1);
	const double& M1_i  = Q.get(i,2);
	const double& M2_i  = Q.get(i,3);
	const double& M3_i  = Q.get(i,4);
	const double  u1_i  = M1_i/rho_i;
	const double  u2_i  = M2_i/rho_i;
	const double  u3_i  = M3_i/rho_i;
	const double& N11_i = Q.get(i,5);
	const double& N12_i = Q.get(i,6);
	const double& N13_i = Q.get(i,7);
	const double& N22_i = Q.get(i,8);
	const double& N23_i = Q.get(i,9);
	const double& N33_i = Q.get(i,10);
	const double  P11_i = N11_i - M1_i*u1_i;
	const double  P12_i = N12_i - M1_i*u2_i;
	const double  P13_i = N13_i - M1_i*u3_i;
	const double  P22_i = N22_i - M2_i*u2_i;
	const double  P23_i = N23_i - M2_i*u3_i;
	const double  P33_i = N33_i - M3_i*u3_i;
	
	const double& rho_e = Q.get(i,11);
	const double& M1_e  = Q.get(i,12);
	const double& M2_e  = Q.get(i,13);
	const double& M3_e  = Q.get(i,14);
	const double  u1_e  = M1_e/rho_e;
	const double  u2_e  = M2_e/rho_e;
	const double  u3_e  = M3_e/rho_e;
	const double& N11_e = Q.get(i,15);
	const double& N12_e = Q.get(i,16);
	const double& N13_e = Q.get(i,17);
	const double& N22_e = Q.get(i,18);
	const double& N23_e = Q.get(i,19);
	const double& N33_e = Q.get(i,20);
	const double  P11_e = N11_e - M1_e*u1_e;
	const double  P12_e = N12_e - M1_e*u2_e;
	const double  P13_e = N13_e - M1_e*u3_e;
	const double  P22_e = N22_e - M2_e*u2_e;
	const double  P23_e = N23_e - M2_e*u3_e;
	const double  P33_e = N33_e - M3_e*u3_e;
	
	const double& B1    = Q.get(i,21);
	const double& B2    = Q.get(i,22);
	const double& B3    = Q.get(i,23);
	const double& E1    = Q.get(i,24);
	const double& E2    = Q.get(i,25);
	const double& E3    = Q.get(i,26);
	    
	// Flux function -- ion fluid equations
	flux.set(i,1,  M1_i );
	flux.set(i,2,  N11_i );
	flux.set(i,3,  N12_i );
	flux.set(i,4,  N13_i );
	flux.set(i,5,  u1_i*N11_i + 2.0*u1_i*P11_i ); 
	flux.set(i,6,  u2_i*N11_i + 2.0*u1_i*P12_i );
	flux.set(i,7,  u3_i*N11_i + 2.0*u1_i*P13_i );
	flux.set(i,8,  u1_i*N22_i + 2.0*u2_i*P12_i );
	flux.set(i,9,  u1_i*N23_i + u3_i*P12_i + u2_i*P13_i );
	flux.set(i,10, u1_i*N33_i + 2.0*u3_i*P13_i );
    
	// Flux function -- electron fluid equations
	flux.set(i,11, M1_e );
	flux.set(i,12, N11_e );
	flux.set(i,13, N12_e );
	flux.set(i,14, N13_e );
	flux.set(i,15, u1_e*N11_e + 2.0*u1_e*P11_e ); 
	flux.set(i,16, u2_e*N11_e + 2.0*u1_e*P12_e );
	flux.set(i,17, u3_e*N11_e + 2.0*u1_e*P13_e );
	flux.set(i,18, u1_e*N22_e + 2.0*u2_e*P12_e );
	flux.set(i,19, u1_e*N23_e + u3_e*P12_e + u2_e*P13_e );
	flux.set(i,20, u1_e*N33_e + 2.0*u3_e*P13_e );
	
	// Flux function -- Maxwell equations
	flux.set(i,21, 0.0 );
	flux.set(i,22, -E3 );
	flux.set(i,23,  E2 );
	flux.set(i,24, 0.0 );
	flux.set(i,25,  cs_light_2*B3 );
	flux.set(i,26, -cs_light_2*B2 );
    }

}
