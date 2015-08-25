#include <cmath>
#include "tensors.h"
#include "TwoFluidParams.h"

// This is a user-supplied routine that projects
// Qvals onto the left eigenvectors ofthe flux 
// Jacobian; the result is stored in Wvals
//
void ProjectLeftEig(const dTensor1& Aux_ave, 
		    const dTensor1& Q_ave,
		    const dTensor2& Qvals,
		    dTensor2& Wvals)
{    
  const int meqn = Qvals.getsize(1);
  const int kmax = Qvals.getsize(2)+1;
  const double gamma = twofluidParams.gamma;
  const double mass_ratio = twofluidParams.mass_ratio;
  const double light_speed = twofluidParams.light_speed;  
  
  // -----------------------------
  // For the IONS
  // -----------------------------
  
  // Average states
  double rho    = Q_ave.get(1);
  double u1     = Q_ave.get(2)/rho;
  double u2     = Q_ave.get(3)/rho;
  double u3     = Q_ave.get(4)/rho;
  double energy = Q_ave.get(5);
  double umag2  = (u1*u1 + u2*u2 + u3*u3);
  double press  = (gamma-1.0e0)*(energy-0.5e0*rho*umag2);
  double c      = sqrt(gamma*press/rho);
  double H      = (energy+press)/rho; 
  
  // Project onto left eigenvectors
  for (int k=1; k<=(kmax-1); k++)
    {
      Wvals.set(1,k, ((umag2/2.0-H-u1*c)*Qvals.get(2,k) 
		      + (umag2/2.0*(c-u1)+H*u1)*Qvals.get(1,k) 
		      + c*(Qvals.get(5,k)-u2*Qvals.get(3,k)
			   - u3*Qvals.get(4,k)))/(c*(2.0*H-umag2)) );
      
      Wvals.set(2,k, 2.0*((H-umag2)*Qvals.get(1,k) + u1*Qvals.get(2,k) 
			  + u2*Qvals.get(3,k) + u3*Qvals.get(4,k) 
			  - Qvals.get(5,k))/(2.0*H-umag2) );
      
      Wvals.set(3,k, Qvals.get(3,k)-u2*Qvals.get(1,k) );
      
      Wvals.set(4,k, Qvals.get(4,k)-u3*Qvals.get(1,k) );
      
      Wvals.set(5,k, ((H-umag2/2.0-u1*c)*Qvals.get(2,k) 
		      + (umag2/2.0*(c+u1)-H*u1)*Qvals.get(1,k)
		      + c*(Qvals.get(5,k)-u2*Qvals.get(3,k) 
			   - u3*Qvals.get(4,k)))/(c*(2.0*H-umag2)) );
    }
  
  // -----------------------------
  // For the ELECTRONS
  // -----------------------------
  
  // Average states
  rho    = Q_ave.get(6);
  u1     = Q_ave.get(7)/rho;
  u2     = Q_ave.get(8)/rho;
  u3     = Q_ave.get(9)/rho;
  energy = Q_ave.get(10);
  umag2  = (u1*u1 + u2*u2 + u3*u3);
  press  = (gamma-1.0e0)*(energy-0.5e0*rho*umag2);
  c      = sqrt(gamma*press/rho);
  H      = (energy+press)/rho; 
  
  // Project onto left eigenvectors
  for (int k=1; k<=(kmax-1); k++)
    {
      Wvals.set(6,k, ((umag2/2.0-H-u1*c)*Qvals.get(7,k) 
		      + (umag2/2.0*(c-u1)+H*u1)*Qvals.get(6,k) 
		      + c*(Qvals.get(10,k)-u2*Qvals.get(8,k)
			   - u3*Qvals.get(9,k)))/(c*(2.0*H-umag2)) );
      
      Wvals.set(7,k, 2.0*((H-umag2)*Qvals.get(6,k) + u1*Qvals.get(7,k) 
			  + u2*Qvals.get(8,k) + u3*Qvals.get(9,k) 
			  - Qvals.get(10,k))/(2.0*H-umag2) );
      
      Wvals.set(8,k, Qvals.get(8,k)-u2*Qvals.get(6,k) );
      
      Wvals.set(9,k, Qvals.get(9,k)-u3*Qvals.get(6,k) );
      
      Wvals.set(10,k, ((H-umag2/2.0-u1*c)*Qvals.get(7,k) 
		       + (umag2/2.0*(c+u1)-H*u1)*Qvals.get(6,k)
		       + c*(Qvals.get(10,k)-u2*Qvals.get(8,k) 
			    - u3*Qvals.get(9,k)))/(c*(2.0*H-umag2)) );
    }
  
  // -----------------------------
  // For the ELECTROMAGNETIC FIELD
  // -----------------------------
    
  // Project onto left eigenvectors
  for (int k=1; k<=(kmax-1); k++)
    {
      Wvals.set(11,k, (light_speed*Qvals.get(12,k) + Qvals.get(16,k))/(2.0*light_speed) );
      
      Wvals.set(12,k, (light_speed*Qvals.get(13,k) - Qvals.get(15,k))/(2.0*light_speed) );
      
      Wvals.set(13,k, Qvals.get(11,k) );
      
      Wvals.set(14,k, Qvals.get(14,k) );
      
      Wvals.set(15,k, (light_speed*Qvals.get(12,k) - Qvals.get(16,k))/(2.0*light_speed) );
      
      Wvals.set(16,k, (light_speed*Qvals.get(13,k) + Qvals.get(15,k))/(2.0*light_speed) );
    }
}
