#include <cmath>
#include "tensors.h"
#include "constants.h"

// This is a user-supplied routine that projects
// Wvals onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Qvals
//
void ProjectRightEig(const dTensor1& Aux_ave, 
		     const dTensor1& Q_ave, 
		     const dTensor2& Wvals,
		     dTensor2& Qvals)
{   
  const int meqn = Qvals.getsize(1);
  const int kmax = Qvals.getsize(2)+1;
  
  // Average states
  const double rho = Q_ave.get(1);
  const double u   = Q_ave.get(2)/rho;
  const double u2  =  u*u;
  const double u3  = u2*u;
  const double p   = Q_ave.get(3)-rho*u2;
  const double q   = Q_ave.get(4)-3.0*p*u-rho*u3;
  const double alpha = 3.0/13.0;
  
  double z1 = u + (13.0/6.0)*q/p - osq13*sqrt(p/rho*(33.0 + 2.0*sqrt(165.0)));
  double z2 = u + (13.0/6.0)*q/p - osq13*sqrt(p/rho*(33.0 - 2.0*sqrt(165.0)));
  double z3 = u + (13.0/6.0)*q/p + osq13*sqrt(p/rho*(33.0 - 2.0*sqrt(165.0)));
  double z4 = u + (13.0/6.0)*q/p + osq13*sqrt(p/rho*(33.0 + 2.0*sqrt(165.0)));

  // Project onto right eigenvectors
  for (int k=1; k<=(kmax-1); k++)
    {
      Qvals.set(1,k, Wvals.get(1,k) + Wvals.get(2,k)
                   + Wvals.get(3,k) + Wvals.get(4,k)  );
      
      Qvals.set(2,k, z1*Wvals.get(1,k) + z2*Wvals.get(2,k) 
		   + z3*Wvals.get(3,k) + z4*Wvals.get(4,k)  );
      
      Qvals.set(3,k, pow(z1,2)*Wvals.get(1,k) + pow(z2,2)*Wvals.get(2,k)
                   + pow(z3,2)*Wvals.get(3,k) + pow(z4,2)*Wvals.get(4,k) );
      
      Qvals.set(4,k, pow(z1,3)*Wvals.get(1,k) + pow(z2,3)*Wvals.get(2,k)
                   + pow(z3,3)*Wvals.get(3,k) + pow(z4,3)*Wvals.get(4,k) );
    }
  

}
