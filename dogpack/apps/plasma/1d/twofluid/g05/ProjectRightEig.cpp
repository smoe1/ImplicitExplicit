#include <cmath>
#include "tensors.h"

// This is a user-supplied routine that projects
// Wvals onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Qvals
//
void ProjectRightEig(const dTensor1& Aux_ave, 
		     const dTensor1& Q_ave, 
		     const dTensor2& Wvals,
		     dTensor2& Qvals)
{    
    int m,k;
    int meqn = Qvals.getsize(1);
    int kmax = Qvals.getsize(2)+1;
    double gamma,rho,u1,u2,u3,energy,umag2,press,c,H;
    double cs_light;

    cs_light = Aux_ave.get(4);
  
    // -----------------------------
    // For the IONS
    // -----------------------------

    // Average states
    gamma  = Aux_ave.get(1);
    rho    = Q_ave.get(1);
    u1     = Q_ave.get(2)/rho;
    u2     = Q_ave.get(3)/rho;
    u3     = Q_ave.get(4)/rho;
    energy = Q_ave.get(5);
    umag2  = (u1*u1 + u2*u2 + u3*u3);
    press  = (gamma-1.0e0)*(energy-0.5e0*rho*umag2);
    c      = sqrt(fabs(gamma*press/rho));
    H      = (energy+press)/rho;  
    
    // Project onto right eigenvectors
    for (k=1; k<=(kmax-1); k++)
    {
        Qvals.set(1,k, Wvals.get(1,k) + Wvals.get(2,k) + Wvals.get(5,k)  );

        Qvals.set(2,k, (u1-c)*Wvals.get(1,k) + u1*Wvals.get(2,k) 
		  + (u1+c)*Wvals.get(5,k) );

	Qvals.set(3,k, u2*(Wvals.get(1,k) + Wvals.get(2,k) + Wvals.get(5,k))
		  + Wvals.get(3,k) );
	
	Qvals.set(4,k, u3*(Wvals.get(1,k) + Wvals.get(2,k) + Wvals.get(5,k))
		  + Wvals.get(4,k) );

	Qvals.set(5,k, (H-u1*c)*Wvals.get(1,k) + umag2/2.0*Wvals.get(2,k)
		  + u2*Wvals.get(3,k) + u3*Wvals.get(4,k) 
		  + (H+u1*c)*Wvals.get(5,k) );
    }


    // -----------------------------
    // For the ELECTRONS
    // -----------------------------

    // Average states
    gamma  = Aux_ave.get(1);
    rho    = Q_ave.get(6);
    u1     = Q_ave.get(7)/rho;
    u2     = Q_ave.get(8)/rho;
    u3     = Q_ave.get(9)/rho;
    energy = Q_ave.get(10);
    umag2  = (u1*u1 + u2*u2 + u3*u3);
    press  = (gamma-1.0e0)*(energy-0.5e0*rho*umag2);
    c      = sqrt(fabs(gamma*press/rho));
    H      = (energy+press)/rho;  
    
    // Project onto right eigenvectors
    for (k=1; k<=(kmax-1); k++)
    {
        Qvals.set(6,k, Wvals.get(6,k) + Wvals.get(7,k) + Wvals.get(10,k)  );

        Qvals.set(7,k, (u1-c)*Wvals.get(6,k) + u1*Wvals.get(7,k) 
		  + (u1+c)*Wvals.get(10,k) );

	Qvals.set(8,k, u2*(Wvals.get(6,k) + Wvals.get(7,k) + Wvals.get(10,k))
		  + Wvals.get(8,k) );
	
	Qvals.set(9,k, u3*(Wvals.get(6,k) + Wvals.get(7,k) + Wvals.get(10,k))
		  + Wvals.get(9,k) );

	Qvals.set(10,k, (H-u1*c)*Wvals.get(6,k) + umag2/2.0*Wvals.get(7,k)
		  + u2*Wvals.get(8,k) + u3*Wvals.get(9,k) 
		  + (H+u1*c)*Wvals.get(10,k) );
    }


    // -----------------------------
    // For the ELECTROMAGNETIC FIELD
    // -----------------------------
    
    // Project onto right eigenvectors
    for (k=1; k<=(kmax-1); k++)
    {
        Qvals.set(11,k, Wvals.get(13,k) );

	Qvals.set(12,k, Wvals.get(11,k) + Wvals.get(15,k) );

	Qvals.set(13,k, Wvals.get(12,k) + Wvals.get(16,k) );
	
	Qvals.set(14,k, Wvals.get(14,k) );

	Qvals.set(15,k, cs_light*( Wvals.get(16,k) - Wvals.get(12,k) ) );

	Qvals.set(16,k, cs_light*( Wvals.get(11,k) - Wvals.get(15,k) ) );
    }
}
